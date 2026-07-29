function run_scheduler_benchmark_v5()
%RUN_SCHEDULER_BENCHMARK_V5 V5.2 MEX scheduler comparison through case57.
% Canonical policies: scan, bandit, novelty. The pool is warmed once and kept
% alive; pool startup and untimed warm-ups are excluded from recorded times.
%
% Optional environment overrides:
%   HEBCPF_SCHED_BENCH_CASES   comma-separated case names
%   HEBCPF_SCHED_BENCH_POLICIES comma-separated policy names
%   HEBCPF_SCHED_BENCH_REPEATS positive integer (default 3)
%   HEBCPF_SCHED_BENCH_WORKERS positive integer; existing pools must match
%   HEBCPF_SCHED_BENCH_ROOT    existing/new output directory

release_root = fileparts(mfilename('fullpath'));
solver_dir = fullfile(release_root,'HEBCPF_MEX_v5.2');
cases_to_run = {'case3','case3TS','case4BB0','case4BBc','case4gs', ...
    'case5loop','case5Salam','case5Salam_mod3','case6ww','case7Salam', ...
    'case9','case9Q','case14mod','case14mod2','case30','case33bw', ...
    'case39','case_ieee30','case57mod','case57'};
policies = {'scan','bandit','novelty'};
repeats = 3;

value = strtrim(getenv('HEBCPF_SCHED_BENCH_CASES'));
if ~isempty(value); cases_to_run = split_csv(value); end
value = strtrim(getenv('HEBCPF_SCHED_BENCH_POLICIES'));
if ~isempty(value); policies = split_csv(value); end
value = strtrim(getenv('HEBCPF_SCHED_BENCH_REPEATS'));
if ~isempty(value); repeats = str2double(value); end
if ~isscalar(repeats) || ~isfinite(repeats) || repeats < 1 || repeats ~= floor(repeats)
    error('HEBCPF:InvalidRepeats','HEBCPF_SCHED_BENCH_REPEATS must be a positive integer.');
end

requested_workers = [];
value = strtrim(getenv('HEBCPF_SCHED_BENCH_WORKERS'));
if ~isempty(value); requested_workers = str2double(value); end
if ~isempty(requested_workers) && (~isscalar(requested_workers) || ...
        ~isfinite(requested_workers) || requested_workers < 1 || ...
        requested_workers ~= floor(requested_workers))
    error('HEBCPF:InvalidWorkers', ...
        'HEBCPF_SCHED_BENCH_WORKERS must be a positive integer.');
end

root_override = strtrim(getenv('HEBCPF_SCHED_BENCH_ROOT'));
if isempty(root_override)
    results_root = fullfile(release_root,'scheduler_benchmark_v5.2',datestr(now,'yyyymmdd_HHMMSS'));
else
    results_root = root_override;
end
if ~exist(results_root,'dir'); mkdir(results_root); end
history_dir = fullfile(results_root,'trace_histories');
reference_dir = fullfile(results_root,'references');
if ~exist(history_dir,'dir'); mkdir(history_dir); end
if ~exist(reference_dir,'dir'); mkdir(reference_dir); end
raw_csv = fullfile(results_root,'scheduler_benchmark_raw.csv');
raw_mat = fullfile(results_root,'scheduler_benchmark_raw.mat');

old_dir = pwd;
cleanup_dir = onCleanup(@() cd(old_dir));
cd(solver_dir);
addpath(solver_dir);
for pi = 1:numel(policies)
    policies{pi} = trace_policy_config(policies{pi});
end
if numel(unique(policies)) ~= numel(policies)
    error('HEBCPF:DuplicatePolicies','Policy list contains aliases or duplicates.');
end

pool_clock = tic;
p = gcp('nocreate');
if isempty(p)
    if isempty(requested_workers)
        p = parpool('local');
    else
        p = parpool('local',requested_workers);
    end
elseif ~isempty(requested_workers) && p.NumWorkers ~= requested_workers
    error('HEBCPF:WorkerCountMismatch', ...
        'Existing pool has %d workers; requested benchmark count is %d.', ...
        p.NumWorkers,requested_workers);
end
pool_start_sec = toc(pool_clock);
fprintf('Scheduler benchmark folder: %s\n',results_root);
fprintf('Workers: %d; pool startup %.3f s (excluded)\n',p.NumWorkers,pool_start_sec);

[~,policy_cfg] = trace_policy_config('bandit');
solver_cfg = solver_params();
metadata = benchmark_metadata(release_root,cases_to_run,policies,repeats, ...
    p.NumWorkers,policy_cfg,solver_cfg,pool_start_sec);
    write_json(fullfile(results_root,'benchmark_metadata_v5.2.json'),metadata);

% Untimed warm-up removes first-use JIT/MEX/pool effects from the three repeats.
global HEBCPOLICY %#ok<GVMIS>
for pi = 1:numel(policies)
    HEBCPOLICY = policies{pi};
    evalc('[warm_result,~]=run_merged_case(''case3'',''parfeval'',false);');
    if ~strcmp(warm_result.status,'PASS')
        error('HEBCPF:WarmupFailed','Warm-up failed for %s: %s',policies{pi},warm_result.message);
    end
end

records = repmat(empty_record(),0,1);
if exist(raw_mat,'file')
    loaded = load(raw_mat,'records');
    if isfield(loaded,'records'); records = loaded.records; end
end

for ci = 1:numel(cases_to_run)
    case_name = cases_to_run{ci};
    reference_file = fullfile(reference_dir,[case_name '_reference.mat']);
    reference_solutions = [];
    if exist(reference_file,'file')
        loaded = load(reference_file,'reference_solutions');
        reference_solutions = loaded.reference_solutions;
    end
    for ri = 1:repeats
        order = mod((ri-1) + (0:numel(policies)-1),numel(policies)) + 1;
        for oi = 1:numel(order)
            policy = policies{order(oi)};
            if record_exists(records,case_name,policy,ri)
                fprintf('SKIP completed: %s / %s / repeat %d\n',case_name,policy,ri);
                continue;
            end
            HEBCPOLICY = policy;
            fprintf('RUN %s / %s / repeat %d\n',case_name,policy,ri);
            run_clock = tic;
            solutions = [];
            console = evalc('[result,solutions]=run_merged_case(case_name,''parfeval'',false);');
            harness_wall = toc(run_clock);

            log_relative = sprintf('%s_%s_r%d.log',case_name,policy,ri);
            log_file = fullfile(results_root,log_relative);
            fid = fopen(log_file,'w');
            if fid >= 0; fwrite(fid,console); fclose(fid); end

            set_distance = NaN;
            if strcmp(result.status,'PASS')
                if isempty(reference_solutions)
                    reference_solutions = solutions;
                    save(reference_file,'reference_solutions','-v7');
                    set_distance = 0;
                elseif isequal(size(reference_solutions),size(solutions))
                    set_distance = order_independent_set_distance(reference_solutions,solutions);
                else
                    set_distance = Inf;
                end
                trace_metrics = result.trace_metrics;
                history_relative = fullfile('trace_histories', ...
                    sprintf('%s_%s_r%d.csv',case_name,policy,ri));
                history_file = fullfile(results_root,history_relative);
                history = trace_history_table(trace_metrics);
                writetable(history,history_file);
            else
                trace_metrics = empty_metrics();
                history_relative = '';
            end

            rec = empty_record();
            rec.case_name = case_name;
            rec.policy = policy;
            rec.repeat = ri;
            rec.status = result.status;
            rec.error_id = result.error_id;
            rec.error_message = result.message;
            rec.workers = p.NumWorkers;
            rec.buses = result.buses;
            rec.solutions = result.solutions;
            rec.total_traces = trace_metrics.total_traces;
            rec.trace_to_50pct = trace_metrics.trace_to_50pct;
            rec.trace_to_90pct = trace_metrics.trace_to_90pct;
            rec.trace_to_99pct = trace_metrics.trace_to_99pct;
            rec.wall_to_90pct_sec = trace_metrics.wall_to_90pct;
            rec.solver_wall_sec = result.wall_time;
            rec.harness_wall_sec = harness_wall;
            rec.discovery_traces = trace_metrics.discovery_traces;
            rec.mean_worker_trace_sec = trace_metrics.mean_worker_trace_sec;
            rec.p90_worker_trace_sec = trace_metrics.p90_worker_trace_sec;
            rec.selection_wall_sec = trace_metrics.selection_wall_sec;
            rec.cursor_checks = trace_metrics.cursor_checks;
            rec.novelty_draws = trace_metrics.novelty_draws;
            rec.max_residual = result.max_residual;
            rec.max_abs_vs_reference = set_distance;
            rec.history_file = strrep(history_relative,'\','/');
            rec.log_file = strrep(log_relative,'\','/');
            records(end+1) = rec; %#ok<AGROW>
            save(raw_mat,'records','cases_to_run','policies','repeats','pool_start_sec','-v7');
            writetable(struct2table(records,'AsArray',true),raw_csv);
            fprintf('DONE %s / %s / r%d: status=%s wall=%.3f solutions=%g traces=%g t90=%g\n', ...
                case_name,policy,ri,result.status,result.wall_time,result.solutions, ...
                trace_metrics.total_traces,trace_metrics.trace_to_90pct);
        end
    end
end

aggregate = aggregate_records(records,cases_to_run,policies);
writetable(aggregate,fullfile(results_root,'scheduler_benchmark_aggregate.csv'));
save(fullfile(results_root,'scheduler_benchmark_complete.mat'), ...
    'records','aggregate','cases_to_run','policies','repeats','pool_start_sec','-v7');
fprintf('Scheduler benchmark complete: %d runs\n',numel(records));
end

function distance = order_independent_set_distance(reference,solutions)
% Symmetric nearest-column L-infinity distance, allowing the HEBC antipode.
distance = 0;
for direction = 1:2
    if direction == 1; from = reference; to = solutions;
    else; from = solutions; to = reference;
    end
    for j = 1:size(from,2)
        delta_pos = max(abs(to-from(:,j)),[],1);
        delta_neg = max(abs(to+from(:,j)),[],1);
        distance = max(distance,min(min(delta_pos,delta_neg)));
    end
end
end

function values = split_csv(value)
values = regexp(value,',','split');
values = cellfun(@strtrim,values,'UniformOutput',false);
values = values(~cellfun(@isempty,values));
end

function yes = record_exists(records,case_name,policy,repeat)
yes = false;
for k = 1:numel(records)
    if strcmp(records(k).case_name,case_name) && strcmp(records(k).policy,policy) && ...
            records(k).repeat == repeat && strcmp(records(k).status,'PASS')
        yes = true;
        return;
    end
end
end

function metrics = empty_metrics()
metrics = struct('total_traces',NaN,'trace_to_50pct',NaN,'trace_to_90pct',NaN, ...
    'trace_to_99pct',NaN,'wall_to_90pct',NaN,'discovery_traces',NaN, ...
    'mean_worker_trace_sec',NaN,'p90_worker_trace_sec',NaN, ...
    'selection_wall_sec',NaN,'cursor_checks',NaN,'novelty_draws',NaN);
end

function rec = empty_record()
rec = struct('case_name','','policy','','repeat',NaN,'status','', ...
    'error_id','','error_message','','workers',NaN,'buses',NaN, ...
    'solutions',NaN,'total_traces',NaN,'trace_to_50pct',NaN, ...
    'trace_to_90pct',NaN,'trace_to_99pct',NaN,'wall_to_90pct_sec',NaN, ...
    'solver_wall_sec',NaN,'harness_wall_sec',NaN,'discovery_traces',NaN, ...
    'mean_worker_trace_sec',NaN,'p90_worker_trace_sec',NaN, ...
    'selection_wall_sec',NaN,'cursor_checks',NaN,'novelty_draws',NaN, ...
    'max_residual',NaN,'max_abs_vs_reference',NaN,'history_file','', ...
    'log_file','');
end

function history = trace_history_table(metrics)
n = numel(metrics.solution_history);
history = table((1:n)',double(metrics.solution_history(:)), ...
    double(metrics.new_solutions_history(:)),metrics.worker_seconds_history(:), ...
    metrics.completion_elapsed_history(:),double(metrics.equation_history(:)), ...
    'VariableNames',{'trace','solutions','new_solutions','worker_seconds', ...
    'completion_elapsed_seconds','equation'});
end

function aggregate = aggregate_records(records,cases_to_run,policies)
rows = repmat(empty_aggregate_row(),0,1);
for ci = 1:numel(cases_to_run)
    for pi = 1:numel(policies)
        idx = arrayfun(@(r) strcmp(r.case_name,cases_to_run{ci}) && ...
            strcmp(r.policy,policies{pi}) && strcmp(r.status,'PASS'),records);
        rr = records(idx);
        if isempty(rr); continue; end
        row = empty_aggregate_row();
        row.case_name = cases_to_run{ci};
        row.policy = policies{pi};
        row.runs = numel(rr);
        row.solutions = mean([rr.solutions]);
        row.total_traces = mean([rr.total_traces]);
        row.trace_to_90pct_mean = mean([rr.trace_to_90pct]);
        row.trace_to_90pct_std = std([rr.trace_to_90pct]);
        row.wall_sec_mean = mean([rr.solver_wall_sec]);
        row.wall_sec_std = std([rr.solver_wall_sec]);
        row.wall_to_90pct_sec_mean = mean([rr.wall_to_90pct_sec]);
        row.mean_worker_trace_sec = mean([rr.mean_worker_trace_sec]);
        row.p90_worker_trace_sec = mean([rr.p90_worker_trace_sec]);
        row.selection_wall_sec_mean = mean([rr.selection_wall_sec]);
        row.discovery_traces_mean = mean([rr.discovery_traces]);
        row.cursor_checks_mean = mean([rr.cursor_checks]);
        row.novelty_draws_mean = mean([rr.novelty_draws]);
        row.max_residual = max([rr.max_residual]);
        row.max_abs_vs_reference = max([rr.max_abs_vs_reference]);
        rows(end+1) = row; %#ok<AGROW>
    end
end
aggregate = struct2table(rows,'AsArray',true);
end

function row = empty_aggregate_row()
row = struct('case_name','','policy','','runs',NaN,'solutions',NaN, ...
    'total_traces',NaN,'trace_to_90pct_mean',NaN,'trace_to_90pct_std',NaN, ...
    'wall_sec_mean',NaN,'wall_sec_std',NaN,'wall_to_90pct_sec_mean',NaN, ...
    'mean_worker_trace_sec',NaN,'p90_worker_trace_sec',NaN, ...
    'selection_wall_sec_mean',NaN,'discovery_traces_mean',NaN, ...
    'cursor_checks_mean',NaN,'novelty_draws_mean',NaN,'max_residual',NaN, ...
    'max_abs_vs_reference',NaN);
end

function metadata = benchmark_metadata(release_root,cases_to_run,policies,repeats, ...
        workers,policy_cfg,solver_cfg,pool_start_sec)
metadata = struct();
metadata.schema_version = 1;
metadata.benchmark_id = datestr(now,'yyyymmdd_HHMMSS');
metadata.suite_version = '5.2.0';
metadata.generated_at = datestr(now,30);
metadata.machine_role = 'benchmark execution machine';
metadata.machine_identifier = strtrim(getenv('COMPUTERNAME'));
if isempty(metadata.machine_identifier); metadata.machine_identifier = 'not recorded'; end
metadata.operating_system = computer;
metadata.matlab_version = version;
metadata.matlab_release = version('-release');
metadata.matpower_version = detect_matpower_version();
metadata.worker_count = workers;
metadata.pool_start_sec_excluded = pool_start_sec;
metadata.repeat_count = repeats;
metadata.cases = cases_to_run;
metadata.policies = policies;
metadata.git_commit = detect_git_commit(release_root);
metadata.timing_convention = ['Pool startup and one untimed warm-up per policy ', ...
    'are excluded; policy order rotates by repeat.'];
metadata.scheduler_parameters = policy_cfg;
metadata.deduplication = struct('seed',7,'signature_count_max',7, ...
    'tolerance',solver_cfg.dedup_tol,'norm','L1','antipodal',true);
metadata.slope_max = solver_cfg.slope_max;
metadata.publication_status = 'local benchmark output; not official until reviewed';
end

function value = detect_matpower_version()
value = 'not detected';
try
    info = mpver;
    if isstruct(info)
        if isfield(info,'Version'); value = char(string(info.Version));
        elseif isfield(info,'version'); value = char(string(info.version));
        end
    else
        value = char(string(info));
    end
catch
end
end

function value = detect_git_commit(release_root)
value = 'not recorded';
command = sprintf('git -C "%s" rev-parse HEAD',release_root);
[status,text] = system(command);
if status == 0 && ~isempty(strtrim(text)); value = strtrim(text); end
end

function write_json(file_name,value)
fid = fopen(file_name,'w');
if fid < 0; error('HEBCPF:MetadataWrite','Cannot write %s.',file_name); end
cleanup = onCleanup(@() fclose(fid));
fwrite(fid,jsonencode(value,'PrettyPrint',true),'char');
fwrite(fid,newline,'char');
end
