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

%% general settings
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

%% R2022a desktop stability fix: do NOT use parallel.pool.Constant.
%  On MATLAB R2022a (Windows, desktop/GUI) the parfeval + parallel.pool.Constant
%  combination can crash MATLAB during heavier cases (e.g. case14mod). The system
%  matrices here are small, so we pass them directly to parfeval (serialised per
%  task) via trace_equation_worker. Verified equivalent; avoids the crash trigger.

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
        futures(eqnN) = parfeval(p, @trace_equation_worker, 3, ...
            MS, M0, Tinv, Ma, Mp, Mq, Ybus, ...
            x_ini, equationnumber_local, solutionnumber, param, ...
            bb, ba, km0, bus_n, I, degree, numofcons, svec_local);
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
            futures(nextSubmit) = parfeval(p, @trace_equation_worker, 3, ...
                MS, M0, Tinv, Ma, Mp, Mq, Ybus, ...
                x_ini, equationnumber_local, solutionnumber, param, ...
                bb, ba, km0, bus_n, I, degree, numofcons, svec_local);
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
            dZoldnew1 = sum(abs(Zsave(variablelist,:) - Z_result(variablelist,eqind1)), 1);
            if Zsave_temp_count > 0
                dZoldnew2 = sum(abs(Zsave_temp(variablelist,1:Zsave_temp_count) - Z_result(variablelist,eqind1)), 1);
                dZoldnew = [dZoldnew1 dZoldnew2];
            else
                dZoldnew = dZoldnew1;
            end
            eqind2 = find(dZoldnew < 4*10^(-7),1);

            ddZoldnew1 = sum(abs(Zsave(variablelist,:) + Z_result(variablelist,eqind1)), 1);
            if Zsave_temp_count > 0
                ddZoldnew2 = sum(abs(Zsave_temp(variablelist,1:Zsave_temp_count) + Z_result(variablelist,eqind1)), 1);
                ddZoldnew = [ddZoldnew1 ddZoldnew2];
            else
                ddZoldnew = ddZoldnew1;
            end
            eqind3 = find(ddZoldnew < 4*10^(-7),1);

            if ~isempty(eqind2)
                New_solu_match{completedIdx}(1,eqind1) = eqind2;
            elseif ~isempty(eqind3)
                New_solu_match{completedIdx}(1,eqind1) = eqind3;
            else
                numberofsolutions = numberofsolutions + 1;
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
    Zsave = [Zsave Zsave_temp(:, 1:Zsave_temp_count)];
    clear futures Zsave_temp;
    solutionnumber = solutionnumber + 1;

    %% save temporary data for a certain period of time
    if totalTime/30000>totalTime_ind
        totalTime_ind = totalTime_ind+1;
%         save case118_temp.mat
    end
end
[Zsave,VBook]=canonicalize_solutions(Zsave,VBook);
clear futures Zsave_temp;
fprintf('Overall executing time: %g\n', toc(tcpu_start));
