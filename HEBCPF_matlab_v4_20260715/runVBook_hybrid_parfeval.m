%% runVBook_hybrid_parfeval.m  (pure MATLAB v4 queue scheduler)
% Global work-queue parfeval driver with per-equation serialization.
% At most one trace for each equation is in flight, preserving VBook pruning
% while avoiding the row barriers used by v3.
%
% RESUME: when the workspace holds a VBook/Zsave checkpoint, the search
% continues from its uncovered pairs instead of re-seeding. Checkpoints are
% written periodically to temp_result.mat.

Tholo=0; Tresol=0; Tcont=0; Tdata=0;
if ~exist('NY','var'); NY=0; NYh=0; end
clear Ysave;

tcpu_start = tic;
numberofequations = length(bb);
equationnumber=1;
N.PV = find(mpc.bus(:,2)==2);
N.PQ = find(mpc.bus(:,2)==1);
N.slack = find(mpc.bus(:,2)==3);
variablelist = 1:numofvar;
numberofvariables = length(variablelist);
sp_ = solver_params(); dedup_tol = sp_.dedup_tol; %#ok<NASGU>
Fail=[];
count = 1;

%% Resume from a compatible checkpoint; otherwise seed the solution set.
initbypass = 0;
if exist('VBook','var')
    VBook = uint32(VBook);   % native class (old double checkpoints converted)
    numberofsolutions = size(VBook,1) %#ok<NOPTS>
    if size(Zsave,2) > numberofsolutions
        Zsave = Zsave(:,1:numberofsolutions);
        initbypass = 1;
        count = max(max(VBook))+1 %#ok<NOPTS>
    elseif size(Zsave,2)==numberofsolutions
        initbypass = 1;
        count = max(max(VBook))+1 %#ok<NOPTS>
    else
        fprintf('\nThe number of saved solutions is less than the bookkeeping count. Start from beginning. \n');
    end
end

if initbypass==1
    if ~exist('totalTime','var');     totalTime = 0;     end
    if ~exist('totalTime_ind','var'); totalTime_ind = 1; end
else
    clear VBook Zsave Vsave;
    Zsave=[]; Vsave=[];
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
    % uint32 from birth: trace ids are small integers; halves VBook RAM
    VBook = zeros(numberofsolutions, 1, 'uint32');
    VBook(1:numberofsolutions,1) = count;
    Zsave = [Zsave Z];
    fprintf('No. Solu: %d, Solu No.: %d, No. Eqt: %d, Eqt No.: %d, Time: %g \n', ...
        numberofsolutions, solutionnumber, numberofequations, 1, Time);
end

%% Deterministic keyed deduplication, matching the v4 MEX driver.
V_dd = size(Zsave,1);
rs_dd = RandStream('mt19937ar','Seed',7);
sigvars_dd = sort(randperm(rs_dd, V_dd, min(7,V_dd)));
dd = KeyedDedup(sigvars_dd, 4e-7, false, max(1024, numberofsolutions), [], 'L1');
dd.seed(Zsave(:, 1:numberofsolutions));

param = solver_params();
param.fac = fac;
param.nsolu = 1;
param.neqt = 1;

p = gcp('nocreate');
if isempty(p); p = parpool('local'); end
numWorkers = p.NumWorkers;

%% Per-equation pending lists: VBook zeros are exactly the uncovered pairs.
pend = cell(numberofequations,1);
pendHead = ones(numberofequations,1);
VBw = size(VBook,2);
for e_ = 1:numberofequations
    if e_ <= VBw
        pend{e_} = find(VBook(1:numberofsolutions,e_)==0);
    else
        pend{e_} = (1:numberofsolutions)';
    end
end
eq_busy = false(numberofequations,1);
eq_ptr = 1;
Zcap = size(Zsave,2);
futs = parallel.FevalFuture.empty(0,1);
fut_s = zeros(0,1); fut_e = zeros(0,1);

while true
    %% Fill worker slots fairly with untraced (solution,equation) pairs.
    while numel(futs) < numWorkers
        found = false;
        for scan = 1:numberofequations
            e_ = mod(eq_ptr + scan - 2, numberofequations) + 1;
            if eq_busy(e_); continue; end
            while pendHead(e_) <= numel(pend{e_})
                s_ = pend{e_}(pendHead(e_));
                pendHead(e_) = pendHead(e_) + 1;
                if size(VBook,1) >= s_ && size(VBook,2) >= e_ && VBook(s_,e_) > 0
                    continue;
                end
                found = true;
                break;
            end
            if found; eq_ptr = e_ + 1; break; end
        end
        if ~found; break; end
        svec_local = zeros(numberofequations,1);
        svec_local(e_) = 1;
        nf = numel(futs) + 1;
        futs(nf,1) = parfeval(p, @trace_equation_worker, 3, ...
            MS, M0, Tinv, Ma, Mp, Mq, Ybus, Zsave(:,s_), e_, s_, param, ...
            bb, ba, km0, bus_n, I, degree, numofcons, svec_local); %#ok<AGROW>
        fut_s(nf,1) = s_; fut_e(nf,1) = e_; %#ok<AGROW>
        eq_busy(e_) = true;
    end
    if isempty(futs); break; end

    try
        [idx, Z_result, Time_result, NY_inc] = fetchNext(futs);
    catch fetchErr
        for fi = 1:numel(futs)
            if isvalid(futs(fi)) && strcmp(futs(fi).State,'finished') && ~isempty(futs(fi).Error)
                fprintf('\n[parfeval error] solution %d, equation %d:\n', fut_s(fi), fut_e(fi));
                disp(getReport(futs(fi).Error));
            end
        end
        error('parfeval task failed: %s', fetchErr.message);
    end
    s_done = fut_s(idx); e_done = fut_e(idx);
    futs(idx) = []; fut_s(idx) = []; fut_e(idx) = [];
    eq_busy(e_done) = false;
    NY = NY + NY_inc;
    totalTime = totalTime + Time_result;
    count = count + 1;
    fprintf('No. Solu: %d, Solu No.: %d, No. Eqt: %d, Eqt No.: %d, Time: %g\n', ...
        numberofsolutions, s_done, numberofequations, e_done, Time_result);

    %% Add genuine discoveries and schedule them on every equation.
    match = zeros(1, size(Z_result,2));
    match(1) = s_done;
    for k = 2:size(Z_result,2)
        [midx, isNew] = dd.query_or_add(Z_result(:,k));
        if ~isNew
            match(k) = midx;
        else
            numberofsolutions = numberofsolutions + 1;
            if numberofsolutions > Zcap
                Zgrow = max(4096, ceil(0.5*Zcap));
                Zsave(:, Zcap+Zgrow) = 0;
                Zcap = Zcap + Zgrow;
            end
            Zsave(:,numberofsolutions) = Z_result(:,k);
            match(k) = numberofsolutions;
            for e_ = 1:numberofequations
                pend{e_}(end+1,1) = numberofsolutions; %#ok<AGROW>
            end
        end
    end
    VBook(match,e_done) = count;

    %% Save only serializable solver state; futures and the pool are rebuilt.
    if totalTime/200000 > totalTime_ind
        totalTime_ind = totalTime_ind + 1;
        % Save only the resume state. Trim unused solution capacity and keep
        % the exact preprocessed problem and ellipse state so main.m need not
        % be rerun. The numerical arrays benefit from uncompressed v7.3 I/O.
        Zsave = Zsave(:, 1:numberofsolutions);
        Zcap = size(Zsave,2);
        save('temp_result.mat', 'Zsave','VBook','mpc','Ybus','bus_n', ...
             'numofvar','numofcons','Mp','Mq','Ma','ba','I','degree', ...
             'MS','M0','T','Tinv','bb','km0','fac','solu','solu0','wait', ...
             'solutionnumber','count','numberofsolutions', ...
             'totalTime','totalTime_ind','NY','NYh', '-v7.3', '-nocompression');
        fprintf('[checkpoint] %d solutions, %d traces\n\n\n', numberofsolutions, count);
    end
end

if size(Zsave,2) > numberofsolutions; Zsave = Zsave(:, 1:numberofsolutions); end
[Zsave,VBook] = canonicalize_solutions(Zsave,VBook);
clear futs;
fprintf('Overall executing time: %g\n', toc(tcpu_start));
