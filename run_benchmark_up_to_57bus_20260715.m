function run_benchmark_up_to_57bus_20260715()
%RUN_BENCHMARK_UP_TO_57BUS_20260715
% Complete parfeval benchmark for the HEBCPF 2026.07.15 v4 release.
%
% Scope:
%   - HEBCPF_MEX_v4_20260715
%   - HEBCPF_matlab_v4_20260715
%   - all bundled cases up to and including 57 buses
%
% Timing convention:
%   The parallel pool is started before the timed case loop. Per-case
%   timings therefore exclude pool startup and shutdown.

release_root = fileparts(mfilename('fullpath'));
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
results_root = fullfile(release_root, 'benchmark_results_20260715', timestamp);
if ~exist(results_root, 'dir')
    mkdir(results_root);
end

cases_to_run = { ...
    'case3', ...
    'case3TS', ...
    'case4BB0', ...
    'case4BBc', ...
    'case4gs', ...
    'case5loop', ...
    'case5Salam', ...
    'case5Salam_mod3', ...
    'case6ww', ...
    'case7Salam', ...
    'case9', ...
    'case9Q', ...
    'case14mod', ...
    'case14mod2', ...
    'case30', ...
    'case33bw', ...
    'case39', ...
    'case_ieee30', ...
    'case57mod', ...
    'case57' ...
    };

solvers = struct( ...
    'name', {'MEX_v4_20260715', 'matlab_v4_20260715'}, ...
    'folder', {'HEBCPF_MEX_v4_20260715', 'HEBCPF_matlab_v4_20260715'}, ...
    'resolveFcn', {'Resolve_with_mex', 'Resolve_no_mex'} ...
    );

case_override = strtrim(getenv('HEBCPF_BENCH_CASES'));
if ~isempty(case_override)
    cases_to_run = split_csv_env(case_override);
end

solver_override = strtrim(getenv('HEBCPF_BENCH_SOLVERS'));
if ~isempty(solver_override)
    wanted_solvers = split_csv_env(solver_override);
    keep = false(size(solvers));
    for ii = 1:numel(solvers)
        keep(ii) = any(strcmpi(solvers(ii).name, wanted_solvers)) || ...
            any(strcmpi(solvers(ii).folder, wanted_solvers));
    end
    solvers = solvers(keep);
end

metadata = struct();
metadata.release_root = release_root;
metadata.results_root = results_root;
metadata.timestamp = timestamp;
metadata.cases_to_run = cases_to_run;
metadata.solvers = solvers;
metadata.matlab_version = version;
metadata.computer = computer;
metadata.timing_note = 'Per-case timings exclude parallel pool startup/shutdown.';

fprintf('HEBCPF 2026.07.15 benchmark: all cases up to 57 buses\n');
fprintf('Results folder: %s\n', results_root);

pool_start = tic;
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local');
end
metadata.pool_start_sec = toc(pool_start);
metadata.pool_num_workers = pool.NumWorkers;
fprintf('Parallel pool ready: %d workers; startup %.3f s (excluded from case timings).\n', ...
    pool.NumWorkers, metadata.pool_start_sec);

all_records = struct([]);

for si = 1:numel(solvers)
    solver = solvers(si);
    solver_dir = fullfile(release_root, solver.folder);
    solver_results_dir = fullfile(results_root, solver.name);
    if ~exist(solver_results_dir, 'dir')
        mkdir(solver_results_dir);
    end

    fprintf('\n============================================================\n');
    fprintf('Solver: %s\n', solver.name);
    fprintf('============================================================\n');

    solver_records = struct([]);

    for ci = 1:numel(cases_to_run)
        case_name = cases_to_run{ci};
        fprintf('\n[%s] Case %d/%d: %s\n', solver.name, ci, numel(cases_to_run), case_name);

        log_file = fullfile(solver_results_dir, [case_name '_console.log']);
        snapshot_file = fullfile(solver_results_dir, [case_name '_snapshot.mat']);

        record = empty_record(solver.name, solver.folder, case_name);
        record.pool_num_workers = pool.NumWorkers;
        record.log_file = log_file;
        record.snapshot_file = snapshot_file;

        old_dir = pwd;
        cd(solver_dir);
        diary(log_file);
        diary on;

        try
            if exist('temp_result.mat', 'file')
                delete('temp_result.mat');
            end

            clear VBook Zsave Vsave Z Y Trace Fail dd futs fut_s fut_e
            clear totalTime totalTime_ind NY NYh Tholo Tresol Tcont Tdata
            clear numberofsolutions numberofequations equationnumber solutionnumber
            clear MS M bb M0 km0 bm0 Mss T Tinv Ma Mp Mq Ybus ba x fac degree
            clear MakeJacobianD quad_form_vals trace_equation_worker
            pctRunOnAll clear MakeJacobianD quad_form_vals trace_equation_worker

            t_case = tic;

            t_preprocess = tic;
            fprintf('Obtain the basic data...\n');
            mpc = feval(case_name);
            factor = 1 + 0.19*(0);
            mpc.bus(:,3:4) = mpc.bus(:,3:4) * factor;
            mpc.gen(:,2:3) = mpc.gen(:,2:3) * factor;
            if factor >= 1
                mpc.gen(:,9) = mpc.gen(:,9) * factor;
                mpc.gen(:,4) = mpc.gen(:,9) * 0.6;
                mpc.gen(:,5) = -mpc.gen(:,9) * 0.6;
            end

            wait = 0; %#ok<NASGU>
            [Ybus, ~, ~] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch);
            [bus_n, ~] = size(mpc.bus);
            [branch_n, ~] = size(mpc.branch);
            [gen_n, ~] = size(mpc.gen);
            numofvar = 2*bus_n - 1;
            numofcons = 2*bus_n - 1;
            maxtrial = 1e4;
            record.preprocess_sec = toc(t_preprocess);
            record.bus_n = bus_n;
            record.branch_n = branch_n;
            record.gen_n = gen_n;
            record.numofvar = numofvar;
            record.numofcons = numofcons;

            t_pf = tic;
            fprintf('Obtain MATPOWER solution...\n');
            result_mpc = runpf(mpc);
            record.matpower_sec = toc(t_pf);
            record.matpower_success = result_mpc.success;
            if result_mpc.success ~= 1
                error('MATPOWER solve failed; could not obtain an initial point.');
            end

            t_quad = tic;
            fprintf('Build quadratic-form matrices...\n');
            [Mp, Mq] = get_quadr_mtrx(Ybus, bus_n);
            [Ma, ba, npg, npd, nq, nv, I] = quadr_matrix(Mp, Mq, mpc, bus_n, gen_n, numofvar); %#ok<ASGLU>
            [solu, solu0] = Solu(result_mpc, I.slack, bus_n); %#ok<ASGLU>
            Err = zeros(numofcons, 1);
            for ii = 1:numofcons
                Err(ii,1) = solu' * Ma{ii} * solu - ba(ii);
            end
            record.quadratic_sec = toc(t_quad);
            record.initial_residual_inf = norm(Err, inf);
            fprintf('Max equation error: %.16g\n', record.initial_residual_inf);

            t_ellipse = tic;
            fprintf('Generate high-dimensional ellipse...\n');
            [MS, M, bb, M0, km0, bm0, succ, Mss, T] = ... %#ok<ASGLU>
                MakeEllipse(Ma, ba, maxtrial, bus_n, gen_n, numofvar, numofcons);
            record.ellipse_sec = toc(t_ellipse);
            record.ellipse_success = succ;
            if succ ~= 1
                error('MakeEllipse failed; could not generate the base ellipse.');
            end

            t_resolve = tic;
            fprintf('Resolve initial point...\n');
            resolve_handle = str2func(solver.resolveFcn);
            [x, flag] = resolve_handle(MS, bb, M0, km0, solu, numofcons, 1); %#ok<NASGU>
            record.resolve_sec = toc(t_resolve);
            record.resolve_flag = flag;

            t_tune = tic;
            fprintf('Tune scaling factor...\n');
            fac = Tune_Fac(MS, M0, km0, x, numofcons); %#ok<NASGU>
            record.tune_sec = toc(t_tune);
            record.fac = fac;

            fprintf('Generate holomorphic parameters...\n');
            Tinv = T \ eye(numofvar); %#ok<NASGU>
            degree.max = 41;
            degree.act = 15;
            degree.denominator = fix(degree.act / 2);
            degree.numerator = degree.act - degree.denominator;

            fprintf('Start parfeval tracing...\n');
            t_trace = tic;
            run('runVBook_hybrid_parfeval.m');
            record.trace_wall_sec = toc(t_trace);

            record.total_case_wall_sec = toc(t_case);
            record.status = 'OK';
            record.error_msg = '';
            record.n_solutions = numberofsolutions;
            record.n_traces_recorded = nnz(VBook);
            record.vbook_rows = size(VBook,1);
            record.vbook_cols = size(VBook,2);
            record.vbook_uncovered = nnz(VBook == 0);
            record.worker_trace_time_sec = local_get_if_exists('totalTime', NaN);
            record.NY = local_get_if_exists('NY', NaN);
            params_for_record = solver_params();
            record.slope_max = params_for_record.slope_max;
            record.dedup_tol = params_for_record.dedup_tol;
            record.solution_abs_sum = sum(abs(Zsave(:)));
            record.solution_real_sum = sum(real(Zsave(:)));
            record.solution_imag_sum = sum(imag(Zsave(:)));
            record.solution_inf_norm = norm(Zsave(:), inf);
            record.mat_snapshot_saved = true;

            save(snapshot_file, '-v7.3', ...
                'record', 'Zsave', 'VBook', 'case_name', 'solver', ...
                'metadata', 'params_for_record');

            fprintf('Benchmark OK: %s / %s | solutions=%d | trace wall=%.4f s | total wall=%.4f s\n', ...
                solver.name, case_name, record.n_solutions, ...
                record.trace_wall_sec, record.total_case_wall_sec);

        catch ME
            record.status = 'ERROR';
            record.error_msg = ME.message;
            record.total_case_wall_sec = local_elapsed_if_started('t_case');
            record.trace_wall_sec = NaN;
            record.mat_snapshot_saved = false;
            fprintf('Benchmark ERROR: %s / %s | %s\n', solver.name, case_name, ME.message);
            try
                save(snapshot_file, '-v7.3', 'record', 'case_name', 'solver', 'metadata');
                record.mat_snapshot_saved = true;
            catch saveErr
                fprintf('Could not save error snapshot: %s\n', saveErr.message);
            end
        end

        diary off;
        cd(old_dir);

        solver_records = append_record(solver_records, record);
        all_records = append_record(all_records, record);
        write_records_csv(solver_records, fullfile(solver_results_dir, 'benchmark_summary.csv'));
        write_records_csv(all_records, fullfile(results_root, 'benchmark_summary_all.csv'));
        save(fullfile(results_root, 'benchmark_records.mat'), '-v7.3', 'metadata', 'all_records');
    end
end

metadata.finished_at = datestr(now, 'yyyy-mm-dd HH:MM:SS');
save(fullfile(results_root, 'benchmark_records.mat'), '-v7.3', 'metadata', 'all_records');
write_records_csv(all_records, fullfile(results_root, 'benchmark_summary_all.csv'));

fprintf('\nBenchmark complete.\n');
fprintf('Combined CSV: %s\n', fullfile(results_root, 'benchmark_summary_all.csv'));
fprintf('Combined MAT: %s\n', fullfile(results_root, 'benchmark_records.mat'));

end

function record = empty_record(solver_name, solver_folder, case_name)
record = struct();
record.solver = solver_name;
record.solver_folder = solver_folder;
record.case_name = case_name;
record.status = 'NOT_RUN';
record.error_msg = '';
record.pool_num_workers = NaN;
record.bus_n = NaN;
record.branch_n = NaN;
record.gen_n = NaN;
record.numofvar = NaN;
record.numofcons = NaN;
record.preprocess_sec = NaN;
record.matpower_sec = NaN;
record.matpower_success = NaN;
record.quadratic_sec = NaN;
record.initial_residual_inf = NaN;
record.ellipse_sec = NaN;
record.ellipse_success = NaN;
record.resolve_sec = NaN;
record.resolve_flag = NaN;
record.tune_sec = NaN;
record.fac = NaN;
record.trace_wall_sec = NaN;
record.total_case_wall_sec = NaN;
record.worker_trace_time_sec = NaN;
record.n_solutions = NaN;
record.n_traces_recorded = NaN;
record.vbook_rows = NaN;
record.vbook_cols = NaN;
record.vbook_uncovered = NaN;
record.NY = NaN;
record.slope_max = NaN;
record.dedup_tol = NaN;
record.solution_abs_sum = NaN;
record.solution_real_sum = NaN;
record.solution_imag_sum = NaN;
record.solution_inf_norm = NaN;
record.mat_snapshot_saved = false;
record.log_file = '';
record.snapshot_file = '';
end

function records = append_record(records, record)
if isempty(records)
    records = record;
else
    records(end+1) = record; %#ok<AGROW>
end
end

function write_records_csv(records, csv_file)
if isempty(records)
    return;
end
T = struct2table(records, 'AsArray', true);
writetable(T, csv_file);
end

function value = local_get_if_exists(var_name, default_value)
if evalin('caller', sprintf('exist(''%s'', ''var'')', var_name))
    value = evalin('caller', var_name);
else
    value = default_value;
end
end

function value = local_elapsed_if_started(var_name)
if evalin('caller', sprintf('exist(''%s'', ''var'')', var_name))
    value = evalin('caller', sprintf('toc(%s)', var_name));
else
    value = NaN;
end
end

function values = split_csv_env(raw)
parts = regexp(raw, ',', 'split');
parts = cellfun(@strtrim, parts, 'UniformOutput', false);
parts = parts(~cellfun(@isempty, parts));
values = parts(:).';
end
