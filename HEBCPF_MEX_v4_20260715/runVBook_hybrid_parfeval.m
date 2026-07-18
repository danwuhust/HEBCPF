%% runVBook_hybrid_parfeval.m  (v4: queue scheduler)
% Global work-queue parfeval driver with per-equation serialization.
%
% Replaces the v3 row-barrier driver (kept in this folder as
% runVBook_hybrid_parfeval_barrier.m): that driver finishes ALL equations of
% solution n before starting solution n+1, so each row ends with a straggler
% tail (measured ~55-61% worker utilization on the large cases at 23
% workers).
%
% This driver keeps one pending list PER EQUATION and enforces at most ONE
% in-flight trace per equation. When popping work for equation e, every
% completed e-trace has already stamped VBook, so the pruning information is
% exactly as complete as the barrier's guarantee: zero redundant traces by
% construction, with no barrier. A rotating scan pointer spreads submissions
% fairly across equations. Parallelism is capped at numberofequations
% (= 2*bus_n-1 >= typical worker counts for the systems of interest).
%
% Measured (23 workers, vs the barrier driver on identical kernels):
% 1.10-1.33x wall, biggest on solution-dense cases; advantage shrinks toward
% parity at small pools (<=8 workers) where the barrier idles less.
%
% The solution SET is unchanged: dedup (KeyedDedup) and
% canonicalize_solutions make the final Zsave identical to the barrier
% driver's.
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
sp_ = solver_params(); dedup_tol = sp_.dedup_tol;
Fail=[];
count = 1;                    % number of traces run

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
    %% first trace on the client (seeds the solution set) -- identical to baseline
    clear VBook Zsave Vsave;
    Zsave=[]; Vsave=[];
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
    % uint32 from birth: trace ids are small integers; halves VBook RAM
    VBook = zeros(numberofsolutions, 1, 'uint32');
    VBook(1:numberofsolutions,1) = count;
    Zsave = [Zsave Z];
    fprintf('No. Solu: %d, Solu No.: %d, No. Eqt: %d, Eqt No.: %d, Time: %g \n', ...
        numberofsolutions, solutionnumber, numberofequations, 1, Time);
end

%% KeyedDedup -- same construction as baseline
V_dd = size(Zsave,1);
rs_dd = RandStream('mt19937ar','Seed',7);
sigvars_dd = sort(randperm(rs_dd, V_dd, min(7,V_dd)));
dd = KeyedDedup(sigvars_dd, 4*10^(-7), false, max(1024, numberofsolutions), [], 'L1');
dd.seed(Zsave(:, 1:numberofsolutions));

param = solver_params();
param.fac   = fac;
param.nsolu = 1;
param.neqt  = 1;

%% pool
p = gcp('nocreate');
if isempty(p)
    p = parpool('local');
end
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
eq_busy = false(numberofequations,1);    % at most one in-flight trace per equation
eq_ptr = 1;                              % rotating fairness pointer

Zcap = size(Zsave,2);                    % headroom-managed capacity of Zsave

futs  = parallel.FevalFuture.empty(0,1);
fut_s = zeros(0,1); fut_e = zeros(0,1);

while true
    %% fill free worker slots: next uncovered pending pair on a non-busy equation
    while numel(futs) < numWorkers
        found = false;
        for scan = 1:numberofequations
            e_ = mod(eq_ptr + scan - 2, numberofequations) + 1;
            if eq_busy(e_); continue; end
            while pendHead(e_) <= numel(pend{e_})
                s_ = pend{e_}(pendHead(e_));
                pendHead(e_) = pendHead(e_) + 1;
                if size(VBook,1) >= s_ && size(VBook,2) >= e_ && VBook(s_,e_) > 0
                    continue;            % covered by a completed e_-trace: prune
                end
                found = true;
                break;
            end
            if found
                eq_ptr = e_ + 1;
                break;
            end
        end
        if ~found; break; end
        svec_local = zeros(numberofequations,1);
        svec_local(e_) = 1;
        % NB: index by numel()+1, NOT end+1 -- after deleting the last element
        % the array is 1x0, where end==1: end+1 would gap-fill slot 1 with a
        % default FevalFuture in state 'unavailable' and break fetchNext.
        nf = numel(futs) + 1;
        futs(nf,1) = parfeval(p, @trace_equation_worker, 3, ...
            MS, M0, Tinv, Ma, Mp, Mq, Ybus, ...
            Zsave(:,s_), e_, s_, param, ...
            bb, ba, km0, bus_n, I, degree, numofcons, svec_local); %#ok<AGROW>
        fut_s(nf,1) = s_; fut_e(nf,1) = e_;                        %#ok<AGROW>
        eq_busy(e_) = true;
    end
    if isempty(futs)
        break;                 % no pending work anywhere, nothing in flight -> done
    end

    %% collect one result
    try
        [idx, Z_result, Time_result, NY_inc] = fetchNext(futs);
    catch fetchErr
        fprintf('\n[diag] numel(futs)=%d numWorkers=%d\n', numel(futs), numWorkers);
        for fi = 1:numel(futs)
            fprintf('[diag] fut %d: s=%d e=%d state=%s\n', fi, fut_s(fi), fut_e(fi), futs(fi).State);
            if isvalid(futs(fi)) && strcmp(futs(fi).State,'finished') && ~isempty(futs(fi).Error)
                fprintf('\n[parfeval error] solu %d eqt %d:\n', fut_s(fi), fut_e(fi));
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

    %% dedup this trace's solutions; enqueue genuinely new ones on every equation
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
            Zsave(:, numberofsolutions) = Z_result(:,k);
            match(k) = numberofsolutions;
            for e_ = 1:numberofequations
                pend{e_}(end+1,1) = numberofsolutions; %#ok<AGROW>
            end
        end
    end
    VBook(match, e_done) = count;

    %% periodic checkpoint (same cadence as the v2/v3 drivers, but enabled).
    %  Saves a truncated, resume-consistent state: load the .mat after the
    %  main.m preprocessing and rerun this script to continue the search.
    if totalTime/200000 > totalTime_ind
        totalTime_ind = totalTime_ind + 1;
        % Save ONLY what resuming the search needs. Zsave is truncated to the
        % confirmed solutions (drops headroom padding); VBook ids are small
        % integers, so uint32 is lossless and half the size. The problem and
        % ellipse state ride along so resume uses the exact preprocessed
        % system without rerunning main.m. -nocompression: solution data barely compresses
        % and MATLAB's gzip is single-threaded (measured 32x slower).
        Zsave = Zsave(:, 1:numberofsolutions);   % drop headroom padding (regrows on demand)
        Zcap = size(Zsave,2);                    % keep capacity tracker in sync
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
