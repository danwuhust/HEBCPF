%% runVBook_hybrid_parfeval.m  (v5: queue scheduler + trace-scheduling policy)
% Global work-queue parfeval driver with per-equation serialization.
%
% Replaces the v3 row-barrier driver (kept in this folder as
% runVBook_hybrid_parfeval_barrier.m): that driver finishes ALL equations of
% solution n before starting solution n+1, so each row ends with a straggler
% tail (measured ~55-61% worker utilization on the large cases at 23
% workers).
%
% This driver enforces at most ONE in-flight trace per equation and finds an
% uncovered start for an equation with a per-equation index cursor over VBook
% (no per-equation pending lists -- see the candidate-selection state below),
% so its memory is O(neq) rather than O(ns*neq). Every completed e-trace has
% stamped VBook, so pruning is exactly as complete as the barrier's guarantee:
% zero redundant traces by construction, with no barrier. Parallelism is
% capped at numberofequations (= 2*bus_n-1 >= typical worker counts).
%
% TRACE-SCHEDULING POLICY (global HEBCPOLICY, see POLICY_README.md):
%   'scan'    -- legacy rotating-fairness equation order + index-cursor start.
%   'bandit'  (default) -- descending learned per-equation gain, most-productive
%             equations first; index-cursor start.
%   'novelty' -- bandit order + bounded-sample novelty start (farthest-from-
%             covered on the equation). Experimental; favours frontier starts,
%             which can cost more per trace. Legacy 'diverse' maps to bandit.
% The cursor removes the O(ns*neq) auxiliary pending lists; VBook itself remains
% O(ns*neq). Completed validation cases matched solution sets across modes.
%
% Measured (23 workers, vs the barrier driver on identical kernels):
% 1.10-1.33x wall, biggest on solution-dense cases; advantage shrinks toward
% parity at small pools (<=8 workers) where the barrier idles less.
%
% RESUME: when the workspace holds a VBook/Zsave checkpoint, the search
% continues from its uncovered pairs instead of re-seeding. Checkpoints are
% written periodically to temp_result.mat, including the learned bandit gains
% (eq_gain_d) so resumed bandit/novelty runs warm-start their equation ranking.

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

%% Optional parallel.pool.Constant for the large read-only operators, so they
%  are sent to each worker ONCE instead of copied on every parfeval task --
%  a large memory/transfer saving on big systems (case118). Off by default:
%  the R2022a parfeval + pool.Constant combination can crash the MATLAB
%  DESKTOP, so opt in (typically under -batch / headless) with
%  `setenv('SOLVER_USE_POOL_CONST','true')`.
useConst = strcmpi(getenv('SOLVER_USE_POOL_CONST'),'true');
if useConst
    MS_c=parallel.pool.Constant(MS);   M0_c=parallel.pool.Constant(M0);
    Tinv_c=parallel.pool.Constant(Tinv); Ma_c=parallel.pool.Constant(Ma);
    Mp_c=parallel.pool.Constant(Mp);   Mq_c=parallel.pool.Constant(Mq);
    Ybus_c=parallel.pool.Constant(Ybus);
    fprintf('[mem] pool.Constant enabled: operators sent to workers once\n');
end

%% Candidate selection state: a per-equation index CURSOR over solutions.
%  Materialising per-equation pending lists (find(VBook(:,e)==0)) is O(ns*neq)
%  -- ~1.8 GB at case118 scale and append-only, a real memory pressure. But the
%  pending list is just the uncovered solution indices in INCREASING order, so
%  we don't need to store it: a single integer cursor per equation walks
%  1..numberofsolutions, skipping covered entries with an O(1) VBook test. That
%  is O(neq) memory (no 1.8 GB), yet visits solutions in the SAME index order
%  as the pending list. Index order == DISCOVERY order, so traces start from
%  early (hub) solutions in numerically benign regions and stay cheap; uniform-
%  random starts instead hit recently-found frontier solutions near folds and
%  trace ~5x slower (measured on case118: 1.7 s vs 9.3 s mean worker time).
cur = ones(1, numberofequations);        % per-equation cursor into 1..numberofsolutions
eq_busy = false(numberofequations,1);    % at most one in-flight trace per equation
eq_ptr = 1;                              % rotating fairness pointer

Zcap = size(Zsave,2);                    % headroom-managed capacity of Zsave

%% trace-scheduling policy (global HEBCPOLICY):
%   'scan'    rotating-fairness equation order + index-cursor solution choice
%             (the v4 benchmarked order).
%   'bandit'  equation-gain bandit order + index-cursor solution choice.
%             DEFAULT. Gain sorting is O(neq log neq); cursor scanning is
%             amortized across the run. Auxiliary selection memory is O(neq).
%   'novelty' bandit equation order + bounded-sample novelty solution choice:
%             start each trace from the discovered solution FARTHEST (in an
%             8-dim JL projection, antipodal-aware) from those already covered
%             on that equation (Law 3 -- avoid redundant re-traces of
%             overlapping curves). EXPERIMENTAL: novelty favours frontier
%             solutions, whose continuation traces are markedly costlier than
%             the low-index hub starts the cursor uses (~5x on case118), so it
%             can be SLOWER in wall time despite higher yield-per-trace.
% Completed regression cases matched solution sets and trace totals. Any mode
% can structurally resume the common checkpoint state; policy-local cursors and
% projections are rebuilt, while eq_gain_d is reused by bandit/novelty.
global HEBCPOLICY %#ok<GVMIS>
[trace_policy,policy_cfg] = trace_policy_config(HEBCPOLICY);
HEBCPOLICY = trace_policy;
use_bandit  = any(strcmp(trace_policy, {'bandit','novelty'}));
use_novelty = strcmp(trace_policy, 'novelty');
checkpoint_trace_interval = str2double(strtrim(getenv('HEBCPF_CHECKPOINT_TRACE_INTERVAL')));
if ~isfinite(checkpoint_trace_interval) || checkpoint_trace_interval < 1 || ...
        checkpoint_trace_interval ~= floor(checkpoint_trace_interval)
    checkpoint_trace_interval = Inf;
end
stop_after_checkpoint = strcmpi(strtrim(getenv('HEBCPF_STOP_AFTER_CHECKPOINT')),'true');
% eq_gain_d: plain bandit state ([] under scan), warm-started from a checkpoint.
if ~exist('eq_gain_d','var'); eq_gain_d = []; end   % exists (empty) for saving; keep a loaded value
if use_bandit
    if numel(eq_gain_d) == numberofequations
        eq_gain_d = eq_gain_d(:)';
        fprintf('[policy] warm-started bandit gains from checkpoint\n');
    else
        eq_gain_d = policy_cfg.bandit_optimistic_gain*ones(1, numberofequations);
    end
end
% novelty sample caps + fixed JL projection (rebuilt each run, so it
% need not be checkpointed); distances use base-MATLAB algebra, not pdist2.
if use_novelty
    nov_kc = policy_cfg.novelty_candidate_cap;
    nov_kv = policy_cfg.novelty_reference_cap;
    nov_budget = policy_cfg.novelty_draw_budget;
    kproj  = min(policy_cfg.novelty_projection_dim, numofvar);
    rngP   = RandStream('mt19937ar','Seed',policy_cfg.novelty_projection_seed);
    P_proj = randn(rngP, numofvar, kproj)/sqrt(kproj);
    rngN   = RandStream('mt19937ar','Seed',policy_cfg.novelty_sampling_seed);
end
% Once discovery STALLS (no new solution for policy_patience traces), bandit and
% bandit and novelty hand off one-way to the scan cursor for the mechanical tail.
% Never fires while discovery continues (e.g. case118).
traces_since_new = 0;
policy_patience  = policy_cfg.stall_sweeps*numberofequations;
fprintf('\n\n[policy] trace scheduling = %s\n\n', trace_policy);

%% Trace-history and scheduler-overhead instrumentation.
if exist('trace_solution_history','var') && exist('trace_history_count','var')
    trace_history_count = double(trace_history_count);
    if ~exist('trace_history_offset','var'); trace_history_offset = 0; end
else
    trace_history_offset = count - 1;
    trace_history_count = 0;
    trace_solution_history = zeros(4096,1,'uint32');
    trace_new_solutions = zeros(4096,1,'uint32');
    trace_worker_seconds = zeros(4096,1);
    trace_completion_elapsed = zeros(4096,1);
    trace_equation_history = zeros(4096,1,'uint32');
    if ~initbypass
        trace_history_offset = 0;
        trace_history_count = 1;
        trace_solution_history(1) = uint32(numberofsolutions);
        trace_new_solutions(1) = uint32(numberofsolutions);
        trace_worker_seconds(1) = Time;
        trace_completion_elapsed(1) = Time;
        trace_equation_history(1) = uint32(1);
    end
end
trace_history_capacity = numel(trace_solution_history);
scheduler_select_seconds = 0;
scheduler_cursor_checks = 0;
scheduler_novelty_draws = 0;
trace_elapsed_base = 0;
if ~initbypass; trace_elapsed_base = Time; end
trace_completion_clock = tic;

futs  = parallel.FevalFuture.empty(0,1);
fut_s = zeros(0,1); fut_e = zeros(0,1);

while true
    %% fill free worker slots. Equation ORDER: scan = rotating fairness,
    %  bandit/novelty = descending learned gain. Solution CHOICE: scan/bandit =
    %  index cursor (low-index/discovery order = cheap traces); novelty =
    %  bounded-sample novelty (farthest-from-covered start), cursor as fallback.
    while numel(futs) < numWorkers
        selection_clock = tic;
        if use_bandit
            [~, eq_order] = sort(eq_gain_d, 'descend');   % O(neq log neq), neq ~ 2*bus_n
        else
            eq_order = mod((eq_ptr-1) + (0:numberofequations-1), numberofequations) + 1;
        end
        found = false;
        nC = size(VBook,2);          % VBook grows columns lazily; e_>nC => all uncovered
        for ii = 1:numberofequations
            e_ = eq_order(ii);
            if eq_busy(e_); continue; end
            s_ = 0;
            if use_novelty
                %% bounded-sample novelty: draw <=nov_budget random solutions,
                %  split into sampled uncovered candidates / covered references,
                %  pick the candidate FARTHEST (antipodal, JL-projected) from the
                %  covered set. O(nov_budget) reads + O(nc*nv*kproj) algebra; no
                %  O(ns) scan, no big allocations.
                cand = zeros(1, nov_kc); nc_ = 0;
                covr = zeros(1, nov_kv); nv_ = 0;
                nd_  = 0;
                while nd_ < nov_budget && (nc_ < nov_kc || nv_ < nov_kv)
                    ss = randi(rngN,numberofsolutions); nd_ = nd_ + 1;
                    if e_ > nC || VBook(ss, e_) == 0
                        if nc_ < nov_kc; nc_ = nc_ + 1; cand(nc_) = ss; end
                    else
                        if nv_ < nov_kv; nv_ = nv_ + 1; covr(nv_) = ss; end
                    end
                end
                scheduler_novelty_draws = scheduler_novelty_draws + nd_;
                if nc_ == 0
                    % none sampled uncovered (equation near-complete): the cursor
                    % is guaranteed to reach any remaining uncovered solution.
                    ce = cur(e_);
                    while ce <= numberofsolutions
                        scheduler_cursor_checks = scheduler_cursor_checks + 1;
                        if e_ > nC || VBook(ce, e_) == 0; s_ = ce; ce = ce + 1; break; end
                        ce = ce + 1;
                    end
                    cur(e_) = ce;
                elseif nv_ == 0
                    s_ = cand(1);   % nothing covered on e yet: any candidate
                else
                    A   = P_proj' * Zsave(:, cand(1:nc_));      % kproj x nc_
                    Bp  = P_proj' * Zsave(:, covr(1:nv_));      % kproj x nv_
                    Bpm = [Bp, -Bp];                            % kproj x 2nv_ (antipodal)
                    D2  = sum(A.^2,1).' - 2*(A.'*Bpm) + sum(Bpm.^2,1);  % nc_ x 2nv_ (||a-b||^2)
                    [~, iw] = max(min(D2, [], 2));              % farthest candidate
                    s_ = cand(iw);
                end
            else
                % index cursor: next uncovered solution in increasing index order
                % (skip covered via the O(1) VBook test; discovery-order = cheap)
                ce = cur(e_);
                while ce <= numberofsolutions
                    scheduler_cursor_checks = scheduler_cursor_checks + 1;
                    if e_ > nC || VBook(ce, e_) == 0; s_ = ce; ce = ce + 1; break; end
                    ce = ce + 1;
                end
                cur(e_) = ce;
            end
            if s_ > 0
                found = true;
                if ~use_bandit; eq_ptr = e_ + 1; end   % advance rotating pointer (scan only)
                break;
            end
        end
        scheduler_select_seconds = scheduler_select_seconds + toc(selection_clock);
        if ~found; break; end
        svec_local = zeros(numberofequations,1);
        svec_local(e_) = 1;
        % NB: index by numel()+1, NOT end+1 -- after deleting the last element
        % the array is 1x0, where end==1: end+1 would gap-fill slot 1 with a
        % default FevalFuture in state 'unavailable' and break fetchNext.
        nf = numel(futs) + 1;
        if useConst
            futs(nf,1) = parfeval(p, @trace_equation_worker_const, 3, ...
                MS_c, M0_c, Tinv_c, Ma_c, Mp_c, Mq_c, Ybus_c, ...
                Zsave(:,s_), e_, s_, param, ...
                bb, ba, km0, bus_n, I, degree, numofcons, svec_local); %#ok<AGROW>
        else
            futs(nf,1) = parfeval(p, @trace_equation_worker, 3, ...
                MS, M0, Tinv, Ma, Mp, Mq, Ybus, ...
                Zsave(:,s_), e_, s_, param, ...
                bb, ba, km0, bus_n, I, degree, numofcons, svec_local); %#ok<AGROW>
        end
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
    gain_new = 0;                        % new solutions this trace (policy signal)
    for k = 2:size(Z_result,2)
        [midx, isNew] = dd.query_or_add(Z_result(:,k));
        if ~isNew
            match(k) = midx;
        else
            numberofsolutions = numberofsolutions + 1;
            gain_new = gain_new + 1;
            if numberofsolutions > Zcap
                Zgrow = max(4096, ceil(0.25*Zcap));   % trimmed headroom (was 0.5) to cap peak Zsave
                Zsave(:, Zcap+Zgrow) = 0;
                Zcap = Zcap + Zgrow;
            end
            Zsave(:, numberofsolutions) = Z_result(:,k);
            match(k) = numberofsolutions;
            % (no per-equation pending list to grow: a new solution is uncovered
            %  on every equation, and the index cursor reaches it in due course)
        end
    end
    VBook(match, e_done) = count;   % stamp this trace's curve as covered on e_done

    %% policy state update (bandit/novelty): discounted equation-gain bandit --
    %  the learned mean of new solutions per trace on this equation. Drives the
    %  equation ordering above. gain_new is the new-solution count this trace.
    if use_bandit
        eq_gain_d(e_done) = policy_cfg.bandit_discount_old*eq_gain_d(e_done) + ...
            policy_cfg.bandit_discount_new*gain_new;
    end

    trace_history_count = trace_history_count + 1;
    if trace_history_count > trace_history_capacity
        trace_history_capacity = trace_history_capacity + max(4096,ceil(0.25*trace_history_capacity));
        trace_solution_history(trace_history_capacity,1) = uint32(0);
        trace_new_solutions(trace_history_capacity,1) = uint32(0);
        trace_worker_seconds(trace_history_capacity,1) = 0;
        trace_completion_elapsed(trace_history_capacity,1) = 0;
        trace_equation_history(trace_history_capacity,1) = uint32(0);
    end
    trace_solution_history(trace_history_count) = uint32(numberofsolutions);
    trace_new_solutions(trace_history_count) = uint32(gain_new);
    trace_worker_seconds(trace_history_count) = Time_result;
    trace_completion_elapsed(trace_history_count) = trace_elapsed_base + toc(trace_completion_clock);
    trace_equation_history(trace_history_count) = uint32(e_done);

    %% bandit/novelty -> scan handoff once discovery stalls: drop to the cheap
    %  rotating-cursor scan for the mechanical tail (gain sort and novelty no
    %  longer help). One-way; only fires after discovery is done.
    if gain_new > 0; traces_since_new = 0; else; traces_since_new = traces_since_new + 1; end
    if use_bandit && traces_since_new > policy_patience
        use_bandit = false; use_novelty = false; eq_ptr = 1;
        fprintf('\n\n[policy] discovery stalled (%d traces, no new) -> scan for the tail at trace %d, %d solutions\n\n\n', ...
            traces_since_new, count, numberofsolutions);
    end

    %% periodic checkpoint (same cadence as the v2/v3 drivers, but enabled).
    %  Saves a truncated, resume-consistent state: load the .mat after the
    %  main.m preprocessing and rerun this script to continue the search.
    forced_checkpoint = isfinite(checkpoint_trace_interval) && ...
        mod(count,checkpoint_trace_interval) == 0;
    if totalTime/200000 > totalTime_ind || forced_checkpoint
        if totalTime/200000 > totalTime_ind
            totalTime_ind = totalTime_ind + 1;
        end
        % Save ONLY what resuming the search needs. Zsave is truncated to the
        % confirmed solutions (drops headroom padding); VBook ids are small
        % integers, so uint32 is lossless and half the size. The problem and
        % ellipse state ride along so resume uses the exact preprocessed
        % system without rerunning main.m. -nocompression: solution data barely compresses
        % and MATLAB's gzip is single-threaded (measured 32x slower).
        Zsave = Zsave(:, 1:numberofsolutions);   % drop headroom padding (regrows on demand)
        Zcap = size(Zsave,2);                    % keep capacity tracker in sync
        trace_solution_history = trace_solution_history(1:trace_history_count);
        trace_new_solutions = trace_new_solutions(1:trace_history_count);
        trace_worker_seconds = trace_worker_seconds(1:trace_history_count);
        trace_completion_elapsed = trace_completion_elapsed(1:trace_history_count);
        trace_equation_history = trace_equation_history(1:trace_history_count);
        trace_history_capacity = trace_history_count;
        save('temp_result.mat', 'Zsave','VBook','mpc','Ybus','bus_n', ...
             'numofvar','numofcons','Mp','Mq','Ma','ba','I','degree', ...
             'MS','M0','T','Tinv','bb','km0','fac','solu','solu0','wait', ...
             'solutionnumber','count','numberofsolutions', ...
             'totalTime','totalTime_ind','NY','NYh','eq_gain_d','trace_policy', ...
             'trace_history_offset','trace_history_count','trace_solution_history', ...
             'trace_new_solutions','trace_worker_seconds','trace_completion_elapsed', ...
             'trace_equation_history', ...
             '-v7.3', '-nocompression');
        fprintf('[checkpoint] %d solutions, %d traces\n\n\n', numberofsolutions, count);
        if stop_after_checkpoint
            cancel(futs);
            fprintf('[checkpoint] controlled stop requested after trace %d\n',count);
            break;
        end
    end
end

if size(Zsave,2) > numberofsolutions; Zsave = Zsave(:, 1:numberofsolutions); end
trace_solution_history = trace_solution_history(1:trace_history_count);
trace_new_solutions = trace_new_solutions(1:trace_history_count);
trace_worker_seconds = trace_worker_seconds(1:trace_history_count);
trace_completion_elapsed = trace_completion_elapsed(1:trace_history_count);
trace_equation_history = trace_equation_history(1:trace_history_count);
[Zsave,VBook] = canonicalize_solutions(Zsave,VBook);
clear futs;
if useConst; clear MS_c M0_c Tinv_c Ma_c Mp_c Mq_c Ybus_c; end   % release worker-side operators
fprintf('Overall executing time: %g\n', toc(tcpu_start));
