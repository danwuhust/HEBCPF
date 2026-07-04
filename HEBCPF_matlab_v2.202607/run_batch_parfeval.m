%% run_batch_parfeval.m (parfeval parallel version)
% Batch-run all test systems using runVBook_hybrid_parfeval.m
% Usage: matlab -batch "run('run_batch_parfeval.m')"
%       or in the MATLAB command window: run('run_batch_parfeval.m')

%% =========================================================
%  Configuration: edit here to choose which systems to test
%% =========================================================
cases_to_run = {
    'case3',         % 6 solutions in 3.95 sec
    'case4gs',       % 6 solutions in 0.37 sec
    'case4BB0',      % 14 solutions in 0.74 sec
    'case4BBc',      % 12 solutions in 0.41 sec
    'case5loop',     % 10 solutions in 0.68 sec
    'case5Salam',    % 10 solutions in 0.59 sec
    'case6ww',       % 6 solutions in 0.37 sec
    'case7Salam',    % 4 solutions in 0.27 sec
    'case9',         % 8 solutions in 0.62 sec
    'case9Q',        % 8 solutions in 0.56 sec
    'case14mod',     % 30 solutions in 3.66 sec
    'case33bw',      % 16 solutions in 8.58 sec
    'case39',        % 176 solutions in 63.51 sec
%     'case30',        % 472 solutions (large; uncomment to stress-test)
};

%% =========================================================
%  Start the parallel pool (requires Parallel Computing Toolbox)
%% =========================================================
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local'); % start with the default number of cores
    fprintf('Parallel pool started, workers: %d\n', pool.NumWorkers);
else
    fprintf('Parallel pool already running, workers: %d\n', pool.NumWorkers);
end

%% =========================================================
%  Initialize result storage
%% =========================================================
n_cases         = length(cases_to_run);
result_case     = cell(n_cases, 1);
result_nsolu    = zeros(n_cases, 1);
result_time     = zeros(n_cases, 1);
result_status   = cell(n_cases, 1);
result_errmsg   = cell(n_cases, 1);

%% =========================================================
%  Main loop: run system by system
%% =========================================================
for ci = 1:n_cases
    case_name = cases_to_run{ci};
    fprintf('\n%s\n', repmat('=', 1, 60));
    fprintf('[%d/%d] system: %s\n', ci, n_cases, case_name);
    fprintf('%s\n', repmat('=', 1, 60));

    result_case{ci}   = case_name;
    result_status{ci} = 'OK';
    result_errmsg{ci} = '';

    try
        % Invalidate case-specific sparse operators on client and workers.
        clear MakeJacobianD quad_form_vals
        pctRunOnAll clear MakeJacobianD quad_form_vals trace_equation_worker
        %% Clear the workspace, keeping loop-control and result variables
        clearvars -except cases_to_run n_cases ci pool ...
                          result_case result_nsolu result_time ...
                          result_status result_errmsg case_name

        %% --------------------------------------------------
        %  Step 1: dynamically load the system and preprocess (equivalent to main.m)
        %% --------------------------------------------------
        fprintf('Obtain the basic data.... \n');
        tic;

        mpc = feval(case_name);

        factor = 1 + 0.19*(0);
        mpc.bus(:,3:4) = mpc.bus(:,3:4) * factor;
        mpc.gen(:,2:3) = mpc.gen(:,2:3) * factor;
        if factor >= 1
            mpc.gen(:,9) = mpc.gen(:,9) * factor;
            mpc.gen(:,4) = mpc.gen(:,9) * 0.6;
            mpc.gen(:,5) = -mpc.gen(:,9) * 0.6;
        end

        wait = 0;

        [Ybus, ~, ~] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch);
        [bus_n, ~]    = size(mpc.bus);
        [branch_n, ~] = size(mpc.branch);
        [gen_n, ~]    = size(mpc.gen);
        numofvar  = 2*bus_n - 1;
        numofcons = 2*bus_n - 1;
        maxtrial  = 1e4;
        toc;

        %% --------------------------------------------------
        %  Step 2: Matpower initial solution
        %% --------------------------------------------------
        tic;
        fprintf('\nobtain solution...\n');
        result_mpc = runpf(mpc);
        toc;

        if result_mpc.success ~= 1
            error('Matpower solve failed; could not obtain an initial point.');
        end

        %% --------------------------------------------------
        %  Step 3: Build quadratic-form matrices
        %% --------------------------------------------------
        tic;
        fprintf('\nobtain quadratic matrices of power balancing equations...\n');
        [Mp, Mq] = get_quadr_mtrx(Ybus, bus_n);
        toc;

        tic;
        fprintf('\nobtain quadratic matrices for PV, PQ and slack buses...\n');
        [Ma, ba, npg, npd, nq, nv, I] = quadr_matrix(Mp, Mq, mpc, bus_n, gen_n, numofvar);
        [solu, solu0] = Solu(result_mpc, I.slack, bus_n);
        Err = zeros(numofcons, 1);
        for i = 1:numofcons
            Err(i,1) = solu' * Ma{i} * solu - ba(i);
        end
        fprintf('Max equation error (should be ~0): %g\n', norm(Err, inf));
        toc;

        %% --------------------------------------------------
        %  Step 4: Generate high-dimensional ellipse
        %% --------------------------------------------------
        tic;
        fprintf('\ngenerate high dimensional ellipses...\n');
        [MS, M, bb, M0, km0, bm0, succ, Mss, T] = ...
            MakeEllipse(Ma, ba, maxtrial, bus_n, gen_n, numofvar, numofcons);

        if succ ~= 1
            error('MakeEllipse failed; could not generate the base ellipse.');
        end
        toc;

        %% --------------------------------------------------
        %  Step 5: Resolve initial point
        %% --------------------------------------------------
        tic;
        fprintf('\nResolve initial point...\n');
        [x, flag] = Resolve_no_mex(MS, bb, M0, km0, solu, numofcons, 1);
        toc;

        %% --------------------------------------------------
        %  Step 6: Tune Jacobian scaling factor
        %% --------------------------------------------------
        tic;
        fprintf('\nTune scaling factor for Jacobian...\n');
        [fac] = Tune_Fac(MS, M0, km0, x, numofcons);
        toc;

        %% --------------------------------------------------
        %  Step 7: Generate holomorphic parameters
        %% --------------------------------------------------
        fprintf('\nGenerate holomorphic parameters...\n');
        Tinv = T \ eye(numofvar);
        degree.max         = 41;
        degree.act         = 15;
        degree.denominator = fix(degree.act / 2);
        degree.numerator   = degree.act - degree.denominator;

        %% --------------------------------------------------
        %  Step 8: run parfeval parallel power-flow all-solution tracing
        %% --------------------------------------------------
        fprintf('\n--- Starting parfeval parallel power-flow all-solution tracing ---\n');
        t_batch_start = tic;

        run('runVBook_hybrid_parfeval.m');

        %% --------------------------------------------------
        %  Step 9: record results
        %% --------------------------------------------------
        result_nsolu(ci) = numberofsolutions;
        result_time(ci)  = toc(t_batch_start);

        fprintf('\n>>> [done] %s | solutions: %d | CPU time: %.2f s\n', ...
            case_name, numberofsolutions, result_time(ci));

    catch ME
        result_status{ci} = 'ERROR';
        result_errmsg{ci} = ME.message;
        result_nsolu(ci)  = -1;
        result_time(ci)   = -1;
        fprintf('\n>>> [error] %s: %s\n', case_name, ME.message);
    end
end

%% =========================================================
%  Summary printout
%% =========================================================
fprintf('\n\n%s\n', repmat('=', 1, 60));
fprintf('              Summary results (parfeval)\n');
fprintf('%s\n', repmat('=', 1, 60));
fprintf('%-15s %12s %10s %8s\n', 'System', 'CPU time(s)', 'solutions', 'Status');
fprintf('%s\n', repmat('-', 1, 50));

for ci = 1:n_cases
    if strcmp(result_status{ci}, 'OK')
        fprintf('%-15s %12.2f %10d %8s\n', ...
            result_case{ci}, result_time(ci), result_nsolu(ci), result_status{ci});
    else
        fprintf('%-15s %12s %10s %8s  [%s]\n', ...
            result_case{ci}, '-', '-', result_status{ci}, result_errmsg{ci});
    end
end

%% =========================================================
%  Save CSV
%% =========================================================
fid = fopen('results_parfeval_summary.csv', 'w', 'n', 'UTF-8');
fprintf(fid, 'case_name,cpu_time_sec,n_solutions,status,error_msg\n');
for ci = 1:n_cases
    safe_err = strrep(result_errmsg{ci}, ',', ';');
    fprintf(fid, '%s,%.4f,%d,%s,%s\n', ...
        result_case{ci}, result_time(ci), result_nsolu(ci), ...
        result_status{ci}, safe_err);
end
fclose(fid);
cleanup_parallel_artifacts(~usejava('desktop'));

fprintf('\nResults saved to results_parfeval_summary.csv\n');
