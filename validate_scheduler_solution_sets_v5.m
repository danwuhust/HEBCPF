function validation = validate_scheduler_solution_sets_v5()
%VALIDATE_SCHEDULER_SOLUTION_SETS_V5 Order-independent policy-set audit.
% Reuses the scan reference saved by the timing benchmark and performs one
% untimed validation run for bandit and novelty on every case through case57.

release_root = fileparts(mfilename('fullpath'));
solver_dir = fullfile(release_root,'HEBCPF_MEX_v5.2');
benchmark_root = strtrim(getenv('HEBCPF_SCHED_BENCH_ROOT'));
if isempty(benchmark_root)
    error('HEBCPF:MissingBenchmarkRoot', ...
        'Set HEBCPF_SCHED_BENCH_ROOT to the completed benchmark directory.');
end
cases_to_run = {'case3','case3TS','case4BB0','case4BBc','case4gs', ...
    'case5loop','case5Salam','case5Salam_mod3','case6ww','case7Salam', ...
    'case9','case9Q','case14mod','case14mod2','case30','case33bw', ...
    'case39','case_ieee30','case57mod','case57'};
policies = {'bandit','novelty'};
value = strtrim(getenv('HEBCPF_SET_VALIDATION_CASES'));
if ~isempty(value); cases_to_run = split_csv(value); end

old_dir = pwd;
cleanup_dir = onCleanup(@() cd(old_dir)); %#ok<NASGU>
cd(solver_dir); addpath(solver_dir);
p = gcp('nocreate');
if isempty(p); p = parpool('local'); end
global HEBCPOLICY %#ok<GVMIS>
rows = repmat(empty_row(),0,1);
out_csv = fullfile(benchmark_root,'scheduler_solution_set_validation.csv');

for ci = 1:numel(cases_to_run)
    case_name = cases_to_run{ci};
    ref_file = fullfile(benchmark_root,'references',[case_name '_reference.mat']);
    loaded = load(ref_file,'reference_solutions');
    reference = loaded.reference_solutions;
    for pi = 1:numel(policies)
        policy = policies{pi};
        HEBCPOLICY = policy;
        fprintf('VALIDATE %s / %s\n',case_name,policy);
        [result,solutions] = run_merged_case(case_name,'parfeval',false);
        row = empty_row();
        row.case_name = case_name;
        row.policy = policy;
        row.status = result.status;
        row.solutions = result.solutions;
        row.reference_solutions = size(reference,2);
        row.max_residual = result.max_residual;
        if strcmp(result.status,'PASS') && isequal(size(reference),size(solutions))
            row.symmetric_nearest_linf = order_independent_set_distance(reference,solutions);
        else
            row.symmetric_nearest_linf = Inf;
        end
        rows(end+1) = row; %#ok<AGROW>
        writetable(struct2table(rows,'AsArray',true),out_csv);
        fprintf('DONE %s / %s: set distance %.3g\n', ...
            case_name,policy,row.symmetric_nearest_linf);
    end
end
validation = struct2table(rows,'AsArray',true);
assert(all(strcmp(validation.status,'PASS')) && ...
    all(validation.symmetric_nearest_linf <= 4e-7), ...
    'HEBCPF:SetValidation','At least one policy solution set failed validation.');
end

function values = split_csv(value)
values = regexp(value,',','split');
values = cellfun(@strtrim,values,'UniformOutput',false);
values = values(~cellfun(@isempty,values));
end

function row = empty_row()
row = struct('case_name','','policy','','status','','solutions',NaN, ...
    'reference_solutions',NaN,'max_residual',NaN, ...
    'symmetric_nearest_linf',NaN);
end

function distance = order_independent_set_distance(reference,solutions)
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
