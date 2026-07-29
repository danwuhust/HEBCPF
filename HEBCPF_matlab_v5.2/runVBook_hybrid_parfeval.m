%% runVBook_hybrid_parfeval.m  (pure MATLAB v5 queue scheduler)
% Global work-queue parfeval driver with per-equation serialization.
% At most one trace for each equation is in flight, preserving VBook pruning
% while avoiding the row barriers used by v3.
%
% TRACE-SCHEDULING POLICY (global HEBCPOLICY, see POLICY_README.md): the order
% equations are visited to fill worker slots is 'bandit' (default descending
% learned per-equation gain -- a bandit) or 'scan' (legacy rotating fairness;
% so the most-productive curves are traced first and solutions were found
% earlier in measured anytime tests). Pending-solution lookup remains O(1),
% but ranking equations requires an O(neq log neq) sort per dispatch. Completed
% validation cases matched solution sets and trace counts; wall time can vary.
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

%% Optional parallel.pool.Constant for the read-only operators (sent to each
%  worker once instead of per task). Off by default -- the R2022a parfeval +
%  pool.Constant combination can crash the desktop; opt in (headless/-batch)
%  with setenv('SOLVER_USE_POOL_CONST','true').
useConst = strcmpi(getenv('SOLVER_USE_POOL_CONST'),'true');
if useConst
    MS_c=parallel.pool.Constant(MS);   M0_c=parallel.pool.Constant(M0);
    Tinv_c=parallel.pool.Constant(Tinv); Ma_c=parallel.pool.Constant(Ma);
    Mp_c=parallel.pool.Constant(Mp);   Mq_c=parallel.pool.Constant(Mq);
    Ybus_c=parallel.pool.Constant(Ybus);
    fprintf('[mem] pool.Constant enabled: operators sent to workers once\n');
end

%% Candidate selection uses one monotonic solution-index cursor per equation.
% This removes the O(ns*neq) auxiliary pending lists. VBook itself remains
% O(ns*neq), and cursor lookup is amortized rather than worst-case O(1).
cur = ones(1,numberofequations);
eq_busy = false(numberofequations,1);
eq_ptr = 1;
Zcap = size(Zsave,2);

%% Policies: scan, bandit (default), and novelty. Legacy v5 'diverse' maps
% to bandit. Completed regression cases matched solution sets and trace totals.
global HEBCPOLICY %#ok<GVMIS>
[trace_policy,policy_cfg] = trace_policy_config(HEBCPOLICY);
HEBCPOLICY = trace_policy;
use_bandit = any(strcmp(trace_policy,{'bandit','novelty'}));
use_novelty = strcmp(trace_policy,'novelty');
checkpoint_trace_interval = str2double(strtrim(getenv('HEBCPF_CHECKPOINT_TRACE_INTERVAL')));
if ~isfinite(checkpoint_trace_interval) || checkpoint_trace_interval < 1 || ...
        checkpoint_trace_interval ~= floor(checkpoint_trace_interval)
    checkpoint_trace_interval = Inf;
end
stop_after_checkpoint = strcmpi(strtrim(getenv('HEBCPF_STOP_AFTER_CHECKPOINT')),'true');
% eq_gain_d: plain discounted per-equation gain bandit state (mean new
% solutions per trace), kept as a plain variable so it saves/loads in
% temp_result.mat by name. Carried across resume to keep the learned ranking.
if ~exist('eq_gain_d','var'); eq_gain_d = []; end
if use_bandit
    if numel(eq_gain_d) == numberofequations
        eq_gain_d = eq_gain_d(:)';
        fprintf('[policy] warm-started bandit gains from checkpoint\n');
    else
        eq_gain_d = policy_cfg.bandit_optimistic_gain*ones(1,numberofequations);
    end
end
if use_novelty
    nov_kc = policy_cfg.novelty_candidate_cap;
    nov_kv = policy_cfg.novelty_reference_cap;
    nov_budget = policy_cfg.novelty_draw_budget;
    kproj = min(policy_cfg.novelty_projection_dim,numofvar);
    rngP = RandStream('mt19937ar','Seed',policy_cfg.novelty_projection_seed);
    P_proj = randn(rngP,numofvar,kproj)/sqrt(kproj);
    rngN = RandStream('mt19937ar','Seed',policy_cfg.novelty_sampling_seed);
end
traces_since_new = 0;
policy_patience = policy_cfg.stall_sweeps*numberofequations;
fprintf('[policy] trace scheduling = %s\n', trace_policy);

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

futs = parallel.FevalFuture.empty(0,1);
fut_s = zeros(0,1); fut_e = zeros(0,1);

while true
    %% Fill worker slots. Equation order is rotating scan or descending gain;
    % solution choice is the cursor or bounded novelty sampling.
    while numel(futs) < numWorkers
        selection_clock = tic;
        if use_bandit
            [~, eq_order] = sort(eq_gain_d, 'descend');
        else
            eq_order = mod((eq_ptr-1) + (0:numberofequations-1), numberofequations) + 1;
        end
        found = false;
        nC = size(VBook,2);
        for ii = 1:numberofequations
            e_ = eq_order(ii);
            if eq_busy(e_); continue; end
            s_ = 0;
            if use_novelty
                cand = zeros(1,nov_kc); nc_ = 0;
                covr = zeros(1,nov_kv); nv_ = 0;
                nd_ = 0;
                while nd_ < nov_budget && (nc_ < nov_kc || nv_ < nov_kv)
                    ss = randi(rngN,numberofsolutions); nd_ = nd_ + 1;
                    if e_ > nC || VBook(ss,e_) == 0
                        if nc_ < nov_kc; nc_ = nc_ + 1; cand(nc_) = ss; end
                    else
                        if nv_ < nov_kv; nv_ = nv_ + 1; covr(nv_) = ss; end
                    end
                end
                scheduler_novelty_draws = scheduler_novelty_draws + nd_;
                if nc_ == 0
                    ce = cur(e_);
                    while ce <= numberofsolutions
                        scheduler_cursor_checks = scheduler_cursor_checks + 1;
                        if e_ > nC || VBook(ce,e_) == 0; s_ = ce; ce = ce + 1; break; end
                        ce = ce + 1;
                    end
                    cur(e_) = ce;
                elseif nv_ == 0
                    s_ = cand(1);
                else
                    A = P_proj' * Zsave(:,cand(1:nc_));
                    Bp = P_proj' * Zsave(:,covr(1:nv_));
                    Bpm = [Bp,-Bp];
                    D2 = sum(A.^2,1).' - 2*(A.'*Bpm) + sum(Bpm.^2,1);
                    [~,iw] = max(min(D2,[],2));
                    s_ = cand(iw);
                end
            else
                ce = cur(e_);
                while ce <= numberofsolutions
                    scheduler_cursor_checks = scheduler_cursor_checks + 1;
                    if e_ > nC || VBook(ce,e_) == 0; s_ = ce; ce = ce + 1; break; end
                    ce = ce + 1;
                end
                cur(e_) = ce;
            end
            if s_ > 0
                found = true;
                if ~use_bandit; eq_ptr = e_ + 1; end
                break;
            end
        end
        scheduler_select_seconds = scheduler_select_seconds + toc(selection_clock);
        if ~found; break; end
        svec_local = zeros(numberofequations,1);
        svec_local(e_) = 1;
        nf = numel(futs) + 1;
        if useConst
            futs(nf,1) = parfeval(p, @trace_equation_worker_const, 3, ...
                MS_c, M0_c, Tinv_c, Ma_c, Mp_c, Mq_c, Ybus_c, Zsave(:,s_), e_, s_, param, ...
                bb, ba, km0, bus_n, I, degree, numofcons, svec_local); %#ok<AGROW>
        else
            futs(nf,1) = parfeval(p, @trace_equation_worker, 3, ...
                MS, M0, Tinv, Ma, Mp, Mq, Ybus, Zsave(:,s_), e_, s_, param, ...
                bb, ba, km0, bus_n, I, degree, numofcons, svec_local); %#ok<AGROW>
        end
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
            Zsave(:,numberofsolutions) = Z_result(:,k);
            match(k) = numberofsolutions;
        end
    end
    VBook(match,e_done) = count;

    %% Policy state update for bandit and novelty.
    %  learned mean of new solutions per trace on this equation, driving the
    %  equation ordering above.
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

    %% Gain-ordered modes hand off to scan after discovery stalls.
    if gain_new > 0; traces_since_new = 0; else; traces_since_new = traces_since_new + 1; end
    if use_bandit && traces_since_new > policy_patience
        use_bandit = false; use_novelty = false; eq_ptr = 1;
        fprintf('[policy] discovery stalled (%d traces, no new) -> scan for the tail at trace %d, %d solutions\n', ...
            traces_since_new, count, numberofsolutions);
    end

    %% Save only serializable solver state; futures and the pool are rebuilt.
    forced_checkpoint = isfinite(checkpoint_trace_interval) && ...
        mod(count,checkpoint_trace_interval) == 0;
    if totalTime/200000 > totalTime_ind || forced_checkpoint
        if totalTime/200000 > totalTime_ind
            totalTime_ind = totalTime_ind + 1;
        end
        % Save only the resume state. Trim unused solution capacity and keep
        % the exact preprocessed problem and ellipse state so main.m need not
        % be rerun. The numerical arrays benefit from uncompressed v7.3 I/O.
        Zsave = Zsave(:, 1:numberofsolutions);
        Zcap = size(Zsave,2);
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
