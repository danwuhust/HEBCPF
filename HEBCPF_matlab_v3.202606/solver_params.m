function param = solver_params()
%SOLVER_PARAMS  Single source of truth for every tunable constant in the
% all-solutions power-flow tracer. Returns a struct `param` consumed by
% hybrid_traceloops -> {holomorphic_hybrid_5, branch_trace, Resolve}.
%
% Only behaviour-controlling constants are here. Pure mathematical literals in the
% formulas (e.g. the 2 in 2*M0*x, the /2 in J*x/2) are NOT parameters and stay in
% the code. `degree` (Padé order) is kept separate because it is baked into the MEX
% binaries (coder.Constant) — changing it requires rebuilding them (see build_mex).
%
% Per-curve mutable fields (fac, nsolu, neqt) are set by hybrid_traceloops, not here.

% ===================== CURVE / OUTER CONTROL ==============================
param.outnum        = 500;      % max holomorphic<->continuation alternations per curve (cap)
param.handback      = 0.5;    % hand-back hysteresis (x aal_min): how far PAST a turning
                                %   point the continuation goes before returning to the
                                %   holomorphic step. Too small -> wrong/non-closing branch;
                                %   too large -> overshoot. 0.5 is the validated value (see README).
param.escape_tol    = 1e-6;     % min|ba| below which the curve starts via continuation
param.closure_tol   = 1e-7;     % ||z-x0||inf below which the curve is declared closed
param.aal_limit_init = 30;      % sentinel for "no pole found yet" (|aal_limit|<this = found)
% --- continuation initial-step cap at hand-back (prevents fold overshoot; see README) ---
param.init_cap_k    = 1500;     % cap hand-back initial astep at this * holo_arc / sstep (Inf = off)

% ===================== HOLOMORPHIC (Pade) EMBEDDING =======================
param.holonum       = 50;       % max Pade steps per outer iteration before forcing continuation
param.holo_arc_max  = 0.01;     % max holomorphic arc length per step
param.holo_arc_init_frac = 0.5; % first-step arc = holo_arc_max * this
param.holo_restart_mult  = 50;  % first holo step after outstep>1: arc = this * aal_min
param.holo_cold_mult     = 3;   % cold-start holo arc = this * aal_min
param.holo_grow_thresh   = 1.5; % grow arc only if last arc > this * arc-before-last
param.holo_shrink_div    = 2;   % otherwise halve the arc (divide by this)
param.hal_default_div    = 20;  % fallback arc = holo_arc_max / this when no holo step succeeded
param.holo_arc_grow      = 1.4; % arc_ratio growth factor while mismatch is small
param.holo_arc_shrink    = 0.75;% arc_ratio shrink factor while mismatch too large
param.holo_unbounded     = 10;  % |arc_ratio*holo_arc| above this => unbounded trace, skip
param.holo_switch_ratio  = 0.1; % holo_arc/holo_arcp <= this => give up holomorphic (flagh=1)
% --- Pade pole estimate: consensus across a FIXED set of buses (was: 2 random buses) ---
param.pole_buses    = 2:6;      % buses sampled for the pole estimate (auto-filtered to <=bus_n)
param.pole_consensus = true;    % pick smallest pole agreed on by >=half the buses (else global-min)
param.pole_cl_tol   = 0.10;     % relative magnitude tolerance for clustering "the same pole"
% --- predictor-corrector restart fraction (krt) when switching to continuation ---
param.krt_interval_div   = 5;   % krt = (dist-to-pole / this) / last-arc
param.krt_max            = 0.45;% krt clamp
param.krt_flagh          = 0.05;% krt when no pole estimate (holomorphic stall)
param.krt_resfail        = 0.03;% krt when the holomorphic correction failed
param.im_tol        = 1e-3;     % imag tolerance for accepting a Pade pole as a real pole
param.mismatch_error = 1e-4;    % max holomorphic constraint mismatch accepted per step

% ===================== NUMERICAL CONTINUATION (pseudo-arclength) ===========
param.smax          = 20000;    % max continuation steps per branch_trace call
param.sstep         = 1e-5;     % base arclength step scale
param.maximumstep   = 8000;     % max step-size multiplier (phase I)
param.maxphase1     = param.maximumstep/2;    % max step in phase I
param.maxphase2     = param.maxphase1/10;     % max step in phase II
param.minimumstep   = 1e-2;     % min step-size multiplier
param.targetsteps   = 3;        % target corrector iterations (drives step adaptation)
param.imax          = 40;       % max corrector Newton iterations per continuation step
param.maxskip       = 50;       % consecutive numerical stalls before Fail=1
param.maxsolu       = 100;      % max solutions collected per continuation loop before Fail=2
param.slope_max     = 1e4;      % max secant slope to accept a turning point as passed
                                %   (1e4 validated; 5e4+ re-breaks case118 hard curves, see README)
param.PS            = 300;      % phase-II (pinv) duration in steps
param.interval      = 1;        % stride for the turning-point finite-difference test
param.dev           = 20;       % console print interval (display only)
param.aal_min       = 1e-6;     % min homotopy-parameter increment (unit for handback etc.)
param.condn         = 1e-8;     % condition-number reference (legacy)
% --- continuation step-size / stall control ---
param.stall_iter_mult   = 3;    % stall if corrector iters > this * targetsteps
param.step_adapt_denom  = 40;   % step growth: astep*(1+(targetsteps-icount)/this)
param.max_adapt         = true; % adaptive per-step cap: astep <= maximumstep * min(slope_prev/slope,1), recomputed fresh each step from the 2-point secant-slope change (NOT compounded). See [[stall-root-cause]].
param.backup_iter_frac  = 0.4;  % iters < this*imax => gentle backup, else aggressive
param.backup_div_small  = 5;    % step shrink on gentle backup
param.backup_div_large  = 10;   % step shrink on aggressive backup
param.zcount_fail_div   = 20;   % step shrink when a crossing solve fails
param.phase2_ef_mult    = 5;    % phase-II residual tolerance = this * ef
param.phase2_ex_mult    = 10;   % phase-II step tolerance = this * ex
param.phase2_trigger_mult = 10; % enter phase II when astep <= this * minimumstep
% --- continuation startup (cold start) ---
param.alpha_init        = 1e-8; % initial homotopy parameter on cold start
param.alpha_retry_div   = 3;    % retry with -alpha/this if first step fails
param.first_step_imax   = 50;   % Newton cap for the cold-start first/second steps
param.second_step_minstep_div = 10; % first-step min astep = minimumstep/this

% ===================== CORRECTOR (Resolve) ================================
param.resolve_maxiter = 10;     % max Newton iterations in Resolve
param.resolve_err     = 2e-10;  % Resolve residual tolerance (||g||)
param.resolve_err1    = 1e-7;   % Resolve step tolerance (||dx||)
param.resolve_fac     = 4;      % Resolve Jacobian scaling (distinct from continuation fac)

% ===================== CONVERGENCE TOLERANCES (shared) ====================
param.ef            = 1e-8;     % continuation residual tolerance ||F||
param.ex            = 1e-7;     % continuation step tolerance ||dx||

% ===================== DEDUPLICATION ======================================
param.dedup_tol     = 4e-7;     % duplicate-solution L1 threshold (full variablelist)

% ===================== PER-CURVE MUTABLE STATE (reset each curve) =========
param.skipindex     = 0;        % running stall counter (reset per curve)
param.Dind          = 0;        % phase flag: 0 = phase I, !=0 = phase II (reset per curve)
end
