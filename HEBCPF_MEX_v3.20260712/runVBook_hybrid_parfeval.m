%runbook.  Script to trace and save traces with bookkeeping.
%  parfeval pipeline version: collect new solutions as soon as fetchNext returns a result
%  tracing and collection overlap; the main thread is not idle while waiting

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
    elseif N.PQ(1)<N.slack && N.PQ(2)>N.slack
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
sp_ = solver_params(); dedup_tol = sp_.dedup_tol;   % duplicate-solution threshold (single source)
Fail=[]; % failure recording book
count = 1; % count the number of traces

Trace=cell(1,1); % trace record
initbypass =0; % if it is 0, starting from the very beginning; if it is 1, starting from the last solution

%% if VBook exists, check if VBook matches Zsave, and start from the existing solution
if exist('VBook','var')
    numberofsolutions = size(VBook,1)
    if size(Zsave,2)==numberofsolutions % Vbook record matches Zsave record
        initbypass = 1; % starting from the last solution
        count = max(max(VBook))+1 % count the latest number of traces and proceed to the next one
    end
end

%%
if initbypass==1
%     solutionnumber = 1;
else % initialize, starting from the first solution and the first trace
    clear VBook Zsave Vsave;

    Zsave=[];
    Vsave=[];

    svec = zeros(size(bb));
    svec(1) = 1;

    solutionnumber = 1;
    totalTime_ind = 1;

    tic;
    hybrid_traceloops_4_with_mex;
    Time=toc;
    NY=NY+size(Y,2);
    totalTime = Time;

    numberofsolutions = size(Z,2);

    VBook(1:numberofsolutions,1) = count;  %bookkeeping.

    Zsave = [Zsave Z];
    fprintf('No. Solu: %d, Solu No.: %d, No. Eqt: %d, Eqt No.: %d, Time: %g \n',numberofsolutions, solutionnumber, numberofequations, 1,Time);
    pause(wait);
end

% [KeyedDedup] O(1)-amortized solution dedup (replaces the O(N) linear Zsave-scan below --
% critical when resuming a large checkpoint). MATCH = L1 over 3-7 random signature vars (tol 4e-7, exact +
% antipodal) so counts are identical to the old scan; the bucket KEY is a sum over a few
% predefined random variables, which keeps buckets tiny (the full-vector sum would cluster
% and degrade to O(N)). Seeded from the existing solutions so dd.n tracks numberofsolutions.
V_dd = size(Zsave,1);
rs_dd = RandStream('mt19937ar','Seed',7);
sigvars_dd = sort(randperm(rs_dd, V_dd, min(7,V_dd)));   % 3-7 predefined random signature vars
dd = KeyedDedup(sigvars_dd, 4*10^(-7), false, max(1024, numberofsolutions), [], 'L1');
dd.seed(Zsave(:, 1:numberofsolutions));

%% tunable constants -- single source of truth in solver_params.m
% (the parfeval workers also rebuild this inside hybrid_traceloops; kept here so the
% client-side first-solution trace and the workers use the identical config.)
param = solver_params();
param.fac   = fac;
param.nsolu = solutionnumber;
param.neqt  = equationnumber;

%% Get the parallel pool
p = gcp('nocreate');
if isempty(p)
    p = parpool('local');
end

%% Gated parallel.pool.Constant (R2022a desktop stability + big-case efficiency).
%  Headless/-batch: wrap large read-only matrices as parallel.pool.Constant so
%  they are sent to each worker only ONCE (much faster for large systems).
%  MATLAB DESKTOP/GUI: the R2022a parfeval+pool.Constant combination can crash
%  MATLAB, so there we pass the matrices directly to parfeval instead.
%  DEFAULT: pass matrices directly to parfeval (no pool.Constant). This is stable
%  on R2022a in BOTH the desktop and -batch, and is as fast as pool.Constant for
%  these problem sizes -- case57mod measured 273s direct vs 287-296s WITH Constant
%  (the latter also intermittently heap-corrupts at pool shutdown, 0xc0000374).
%  Opt in to pool.Constant only for very large systems via SOLVER_USE_POOL_CONST=true.
useConst = strcmpi(getenv('SOLVER_USE_POOL_CONST'),'true');
if useConst
    MS_c   = parallel.pool.Constant(MS);
    M0_c   = parallel.pool.Constant(M0);
    Tinv_c = parallel.pool.Constant(Tinv);
    Ma_c   = parallel.pool.Constant(Ma);
    Mp_c   = parallel.pool.Constant(Mp);
    Mq_c   = parallel.pool.Constant(Mq);
    Ybus_c = parallel.pool.Constant(Ybus);
end

while numberofsolutions+1>solutionnumber
    equationqueue = setdiff(1:numberofequations,find(VBook(solutionnumber,:)>0));
    Neqnqueue = size(equationqueue,2); % length of equation queue
    New_solu_match = cell(Neqnqueue,1);
    Zsave_ind = 1:numberofsolutions;

    %% trace different loops - sliding-window parfeval
    %  submit at most numWorkers tasks at a time; as each finishes, collect it and submit another
    %  limit the number of simultaneous futures to reduce memory pressure
    x_ini = Zsave(:,solutionnumber);
    numWorkers = p.NumWorkers;
    batchSize = min(numWorkers, Neqnqueue);

    %% submit the first batch
    futures(1:batchSize) = parallel.FevalFuture;
    for eqnN = 1:batchSize
        equationnumber_local = equationqueue(eqnN);
        svec_local = zeros(numberofequations, 1);
        svec_local(equationnumber_local) = 1;
        if useConst
            futures(eqnN) = parfeval(p, @trace_equation_worker_const, 3, ...
                MS_c, M0_c, Tinv_c, Ma_c, Mp_c, Mq_c, Ybus_c, ...
                x_ini, equationnumber_local, solutionnumber, param, ...
                bb, ba, km0, bus_n, I, degree, numofcons, svec_local);
        else
            futures(eqnN) = parfeval(p, @trace_equation_worker, 3, ...
                MS, M0, Tinv, Ma, Mp, Mq, Ybus, ...
                x_ini, equationnumber_local, solutionnumber, param, ...
                bb, ba, km0, bus_n, I, degree, numofcons, svec_local);
        end
    end
    nextSubmit = batchSize + 1;  % index of the next equation to submit

    %% pipeline: collect new solutions as soon as a result arrives, and submit the next task
    Zsave_temp_chunk = 100;
    Zsave_temp = zeros(size(Zsave,1), Zsave_temp_chunk);
    Zsave_temp_count = 0;
    nCompleted = 0;
    while nCompleted < Neqnqueue
        try
            [completedIdx, Z_result, Time_result, NY_inc] = fetchNext(futures);
        catch fetchErr
            for fi = 1:length(futures)
                if isvalid(futures(fi)) && strcmp(futures(fi).State, 'finished') && ~isempty(futures(fi).Error)
                    fprintf('\n[parfeval error] eqt %d:\n', equationqueue(fi));
                    disp(getReport(futures(fi).Error));
                end
            end
            error('parfeval task failed: %s', fetchErr.message);
        end
        nCompleted = nCompleted + 1;

        %% immediately submit the next task (keep workers fully loaded)
        if nextSubmit <= Neqnqueue
            equationnumber_local = equationqueue(nextSubmit);
            svec_local = zeros(numberofequations, 1);
            svec_local(equationnumber_local) = 1;
            if useConst
                futures(nextSubmit) = parfeval(p, @trace_equation_worker_const, 3, ...
                    MS_c, M0_c, Tinv_c, Ma_c, Mp_c, Mq_c, Ybus_c, ...
                    x_ini, equationnumber_local, solutionnumber, param, ...
                    bb, ba, km0, bus_n, I, degree, numofcons, svec_local);
            else
                futures(nextSubmit) = parfeval(p, @trace_equation_worker, 3, ...
                    MS, M0, Tinv, Ma, Mp, Mq, Ybus, ...
                    x_ini, equationnumber_local, solutionnumber, param, ...
                    bb, ba, km0, bus_n, I, degree, numofcons, svec_local);
            end
            nextSubmit = nextSubmit + 1;
        end

        NY = NY + NY_inc;
        totalTime = totalTime + Time_result;

        equationnumber_done = equationqueue(completedIdx);
        fprintf('No. Solu: %d, Solu No.: %d, No. Eqt: %d, Eqt No.: %d, Time: %g\n\n', ...
            numberofsolutions, solutionnumber, numberofequations, equationnumber_done, Time_result);

        %% immediately collect this equation's new solutions (in parallel with other tracing tasks)
        New_solu_match{completedIdx}(1,1) = solutionnumber;
        for eqind1 = 2:size(Z_result, 2)
            % KeyedDedup: one bucketed lookup against ALL solutions found so far
            % (confirmed + this round's), replacing the O(N) linear Zsave/Zsave_temp scan.
            [idx, isNew] = dd.query_or_add(Z_result(:,eqind1));
            if ~isNew
                New_solu_match{completedIdx}(1,eqind1) = idx;
            else
                numberofsolutions = numberofsolutions + 1;   % == idx == dd.n
                Zsave_temp_count = Zsave_temp_count + 1;
                if Zsave_temp_count > size(Zsave_temp, 2)
                    Zsave_temp(:,2*size(Zsave_temp,2))=0;
                end
                Zsave_temp(:, Zsave_temp_count) = Z_result(:,eqind1);
                New_solu_match{completedIdx}(1,eqind1) = numberofsolutions;
            end
        end
        count = count + 1;
        VBook(New_solu_match{completedIdx}, equationqueue(completedIdx)) = count;
    end
    % Append new solutions in place (headroom buffer, grown geometrically & rarely) -- an
    % O(k) write instead of an O(N) full-array copy each iteration. Truncated before
    % canonicalize/save. Zero-yield iterations skipped.
    if Zsave_temp_count > 0
        if numberofsolutions > size(Zsave,2)
            if ~exist('Zsave_grow','var'); Zsave_grow = max(1024, round(0.01*numberofsolutions)); end
            Zsave(:, max(size(Zsave,2)+Zsave_grow, numberofsolutions)) = 0;
            Zsave_grow = 2*Zsave_grow;
        end
        Zsave(:, numberofsolutions-Zsave_temp_count+1 : numberofsolutions) = Zsave_temp(:, 1:Zsave_temp_count);
    end
    clear futures Zsave_temp;
    solutionnumber = solutionnumber + 1;

    %% save temporary data for a certain period of time
    if totalTime/30000>totalTime_ind
        totalTime_ind = totalTime_ind+1;
%         save case118_temp.mat
    end
end
if size(Zsave,2) > numberofsolutions; Zsave = Zsave(:, 1:numberofsolutions); end
[Zsave,VBook]=canonicalize_solutions(Zsave,VBook);
clear futures Zsave_temp;
if useConst; clear MS_c M0_c Tinv_c Ma_c Mp_c Mq_c Ybus_c; end
fprintf('Overall executing time: %g\n', toc(tcpu_start));
