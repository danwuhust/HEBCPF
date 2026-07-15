%runbook.  Script to trace and save traces with bookkeeping.
%Initialize:  start with solution x.
%  parfor version (function-call mode): uses trace_equation_worker instead of an inline script
%  avoids heap corruption caused by workspace bloat

% global bus_n gen_n numofcons numofvar;
Tholo=0;
Tresol=0;
Tcont=0;
Tdata=0;
% totalTime=0;
if exist('NY','var')

else
    NY=0;
    NYh=0;
end
clear Ysave;

tcpu_start = tic; % recording the starting time (wall-clock)
numberofequations = length(bb);
% solutionnumber = 1;
equationnumber=1;
N.PV = find(mpc.bus(:,2)==2);
N.PQ = find(mpc.bus(:,2)==1);
N.slack = find(mpc.bus(:,2)==3);
%{
if size(N.PQ,1)>=2
    if N.PQ(1)<N.slack && N.PQ(2)<N.slack
        variablelist = [N.PQ(1), N.PQ(2)] + bus_n;
    elseif N.PQ(1)<N.slack && N.PQ(2)>N.slash
        variablelist = [N.PQ(1)+bus_n, N.PQ(2) + bus_n-1];
    elseif N.PQ(1)>N.slack && N.PQ(2)<N.slack
        variablelist = [N.PQ(1)+bus_n-1, N.PQ(2) + bus_n];
    elseif N.PQ(1)>N.slack && N.PQ(2)>N.slack
        variablelist = [N.PQ(1), N.PQ(2)] + bus_n -1;
    end
else
    variablelist = bus_n+1:2*bus_n-1;
end
%}
variablelist = 1:numofvar;
numberofvariables = length(variablelist);
Fail=[]; % failure recording book
count = 1; % count the number of traces

Trace=cell(1,1); % trace record
initbypass =0; % if it is 0, starting from the very beginning; if it is 1, starting from the last solution

%% if VBook exists, check if VBook matches Zsave, and start from the existing solution
if exist('VBook','var')
    numberofsolutions = size(VBook,1)
    if size(Zsave,2) > numberofsolutions 
        Zsave = Zsave(:, 1:numberofsolutions); 
        initbypass = 1;                   % resume from the checkpoint
        count = max(max(VBook))+1 %#ok<NOPTS> % next trace id
    elseif size(Zsave,2)==numberofsolutions   % VBook record matches Zsave record
        initbypass = 1;                   % resume from the checkpoint
        count = max(max(VBook))+1 %#ok<NOPTS> % next trace id
    else
        fprintf('\nThe number of saved solutions is less than the bookkeeping count. Start from beginning. \n');
    end
end

%%
if initbypass==1
%     solutionnumber = 1;
    if ~exist('totalTime','var');     totalTime = 0;     end
    if ~exist('totalTime_ind','var'); totalTime_ind = 1; end
else % initialize, starting from the first solution and the first trace
    clear VBook Zsave Vsave;

    Zsave=[];
    Vsave=[];

    svec = zeros(size(bb));
    svec(1) = 1;

    solutionnumber = 1;
    totalTime_ind = 1;

    tic;
    hybrid_traceloops_4_no_mex;
    Time=toc;
    NY=NY+size(Y,2);
    totalTime = Time;

    numberofsolutions = size(Z,2);

    VBook(1:numberofsolutions,1) = count;  %bookkeeping.

    Zsave = [Zsave Z];
    fprintf('No. Solu: %d, Solu No.: %d, No. Eqt: %d, Eqt No.: %d, Time: %g \n',numberofsolutions, solutionnumber, numberofequations, 1,Time);
    pause(wait);
end

% [KeyedDedup] O(1)-amortized solution dedup (replaces the O(N) linear Zsave-scan below).
% MATCH = L1 over 3-5 random signature variables (tol 4e-7); the sum(sig) key keeps buckets
% tiny. Seeded from the existing solutions so dd.n tracks numberofsolutions.
V_dd = size(Zsave,1);
rs_dd = RandStream('mt19937ar','Seed',7);
sigvars_dd = sort(randperm(rs_dd, V_dd, min(7,V_dd)));   % up to 7 predefined random signature vars
dd = KeyedDedup(sigvars_dd, 4*10^(-7), false, max(1024, numberofsolutions), [], 'L1');
dd.seed(Zsave(:, 1:numberofsolutions));

%% general settings
param.sstep = 1e-5;
param.targetsteps = 3;                  % target number of iterations.
param.imax = 40;                        % maximum number of Newton iterations
param.smax = 20000;                     % maximum number of steps
param.maximumstep = 8000;               % maximum step length
param.maxphase1 = param.maximumstep/2;  % maximum step length for phase I
param.maxphase2 = param.maxphase1/10;   % maximum step length for phase II
param.minimumstep = 1e-2;               % minimum step length
param.dev=20;                            % denominator of interval printing
param.skipindex=0;                      % skip index
param.ef=1e-8;                          % tolerence for function value
param.ex=1e-7;                          % tolerence for dx value
param.Dind=0;                           % trigger for changing phases
param.PS=300;                           % phase II total steps
param.condn=1e-8;                       % condition number
param.maxsolu=100;                      % maximum number of solutions for each path
param.maxskip=50;                       % maximum number of numerical stalling
param.fac=fac;
warning off MATLAB:nearlySingularMatrix;

%% special settings
param.mismatch_error=0.1*10^-3; % deafault 0.1*10^-3
param.holonum=50;               % default 50
param.holo_arc_max=0.01;        % default 0.02
param.outnum=500;               % default 500
param.slope_max=4*10^5;         % v4 release default
param.interval=1;               % default 1
param.aal_min=0.1*10^(-5);      % default 0.5*10^(-5)
param.im_tol=1*10^(-3);         % default 1*10^(-3)
param.nsolu=solutionnumber;
param.neqt=equationnumber;

%%
%% Optional worker-local constants. Disabled by default for desktop stability.
useConst = strcmpi(getenv('SOLVER_USE_POOL_CONST'),'true');
if useConst
    MS_c=parallel.pool.Constant(MS); M0_c=parallel.pool.Constant(M0);
    Tinv_c=parallel.pool.Constant(Tinv); Ma_c=parallel.pool.Constant(Ma);
    Mp_c=parallel.pool.Constant(Mp); Mq_c=parallel.pool.Constant(Mq);
    Ybus_c=parallel.pool.Constant(Ybus); getv=@(c) c.Value;
else
    MS_c=MS; M0_c=M0; Tinv_c=Tinv; Ma_c=Ma; Mp_c=Mp; Mq_c=Mq; Ybus_c=Ybus;
    getv=@(c) c;
end

while numberofsolutions+1>solutionnumber
    equationqueue = setdiff(1:numberofequations,find(VBook(solutionnumber,:)>0));
    Neqnqueue = size(equationqueue,2); % length of equation queue
    New_solu = cell(Neqnqueue,1); % new solution cell
    New_solu_match = cell(Neqnqueue,1);
    Zsave_ind = 1:numberofsolutions;

    %% trace different loops - use a function call instead of an inline script
    %  a function call has its own workspace released on completion, giving better memory isolation
    %  avoids the workspace bloat and heap corruption caused by inlining a 400-line script
    x_ini = Zsave(:,solutionnumber);
    NY_inc_arr = zeros(Neqnqueue, 1);   % collect each equation's NY increment
    Time_arr   = zeros(Neqnqueue, 1);   % collect each equation's elapsed time
    parfor eqnN = 1:Neqnqueue
        equationnumber = equationqueue(eqnN);
        svec_local = zeros(numberofequations, 1);
        svec_local(equationnumber) = 1;

        [Z_result, Time_result, NY_inc_result] = trace_equation_worker( ...
            getv(MS_c), getv(M0_c), getv(Tinv_c), getv(Ma_c), ...
            getv(Mp_c), getv(Mq_c), getv(Ybus_c), ...
            x_ini, equationnumber, solutionnumber, param, ...
            bb, ba, km0, bus_n, I, degree, numofcons, svec_local);

        New_solu{eqnN, 1} = Z_result;
        NY_inc_arr(eqnN) = NY_inc_result;
        Time_arr(eqnN) = Time_result;
    end
    %% Aggregate parfor results
    NY = NY + sum(NY_inc_arr);
    totalTime = totalTime + sum(Time_arr);
    for eqnN = 1:Neqnqueue
        fprintf('No. Solu: %d, Solu No.: %d, No. Eqt: %d, Eqt No.: %d, Time: %g\n\n', ...
            numberofsolutions, solutionnumber, numberofequations, equationqueue(eqnN), Time_arr(eqnN));
    end

    %% collecting new solutions (chunked preallocation to reduce memory fragmentation)
    Zsave_temp_chunk = 100;
    Zsave_temp = zeros(size(Zsave,1), Zsave_temp_chunk);
    Zsave_temp_count = 0;
    for eqind0 = 1:Neqnqueue
        New_solu_match{eqind0}(1,1)=solutionnumber;
        for eqind1 = 2:size(New_solu{eqind0},2)
            % KeyedDedup: one bucketed lookup replacing the O(N) linear Zsave/Zsave_temp scan.
            [idx, isNew] = dd.query_or_add(New_solu{eqind0}(:,eqind1));
            if ~isNew
                New_solu_match{eqind0}(1,eqind1)=idx;
            else
                numberofsolutions = numberofsolutions+1;
                Zsave_temp_count = Zsave_temp_count + 1;
                if Zsave_temp_count > size(Zsave_temp, 2)
                    Zsave_temp(:,2*size(Zsave_temp,2))=0;
                end
                Zsave_temp(:, Zsave_temp_count) = New_solu{eqind0}(:,eqind1);
                New_solu_match{eqind0}(1,eqind1)=numberofsolutions;
            end
        end
        count = count + 1;
        VBook(New_solu_match{eqind0},equationqueue(eqind0))=count;
    end
    % Append new solutions in place (headroom buffer, grown geometrically & rarely) -- an
    % O(k) write instead of an O(N) full-array copy each iteration. Truncated before
    % canonicalize/save. Zero-yield iterations skipped.
    if Zsave_temp_count > 0
        if numberofsolutions > size(Zsave,2)
            if ~exist('Zsave_grow','var'); Zsave_grow = max(4096, round(0.01*numberofsolutions)); end
            Zsave(:, max(size(Zsave,2)+Zsave_grow, numberofsolutions)) = 0;
            Zsave_grow = 2*Zsave_grow;
        end
        Zsave(:, numberofsolutions-Zsave_temp_count+1 : numberofsolutions) = Zsave_temp(:, 1:Zsave_temp_count);
    end
    clear Zsave_temp New_solu;
    solutionnumber = solutionnumber +1;

    %% save temporary data for a certain period of time
    if totalTime/30000>totalTime_ind
        totalTime_ind = totalTime_ind+1;
        save('temp_result.mat', '-v7.3', '-regexp', '^(?!(futs|p)$).');
        fprintf('[checkpoint] %d solutions, %d traces\n\n\n', numberofsolutions, count);
    end
end
if size(Zsave,2) > numberofsolutions; Zsave = Zsave(:, 1:numberofsolutions); end
[Zsave,VBook]=canonicalize_solutions(Zsave,VBook);
fprintf('Overall executing time: %g\n', toc(tcpu_start));
