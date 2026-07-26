function results = smoke_test_v5()
%SMOKE_TEST_V5 Final scheduler/implementation smoke matrix.
% Runs case3, case5loop, and case14mod with all three schedulers in both the
% MEX and pure-MATLAB implementations. A persistent pool excludes pool startup
% from repeated setup and keeps this test short enough for release checks.

release_root = fileparts(mfilename('fullpath'));
implementations = {'HEBCPF_MEX_v5','HEBCPF_matlab_v5'};
cases = {'case3','case5loop','case14mod'};
expected = [6,10,30];
policies = {'scan','bandit','novelty'};
rows = repmat(empty_row(),0,1);

p = gcp('nocreate');
if isempty(p); p = parpool('local'); end
for ii = 1:numel(implementations)
    solver_dir = fullfile(release_root,implementations{ii});
    old_dir = pwd;
    cd(solver_dir);
    addpath(solver_dir);
    for ci = 1:numel(cases)
        for pi = 1:numel(policies)
            global HEBCPOLICY %#ok<GVMIS>
            HEBCPOLICY = policies{pi};
            [result,~] = run_merged_case(cases{ci},'parfeval',false);
            row = empty_row();
            row.implementation = implementations{ii};
            row.case_name = cases{ci};
            row.policy = policies{pi};
            row.solutions = result.solutions;
            row.expected_solutions = expected(ci);
            row.max_residual = result.max_residual;
            row.wall_time_sec = result.wall_time;
            if isfield(result.trace_metrics,'total_traces')
                row.total_traces = result.trace_metrics.total_traces;
                row.trace_to_90pct = result.trace_metrics.trace_to_90pct;
            end
            passed = strcmp(result.status,'PASS') && ...
                strcmp(result.scheduler,policies{pi}) && ...
                result.solutions == expected(ci) && ...
                isfinite(result.max_residual) && result.max_residual <= 4e-7;
            row.status = ternary(passed,'PASS','FAIL');
            row.message = result.message;
            rows(end+1) = row; %#ok<AGROW>
            fprintf('%-21s %-9s %-7s %s: %g solutions, %.3g residual\n', ...
                implementations{ii},cases{ci},policies{pi},row.status, ...
                row.solutions,row.max_residual);
        end
    end
    rmpath(solver_dir);
    cd(old_dir);
    clear run_merged_case trace_policy_config
end

results = struct2table(rows,'AsArray',true);
out_dir = fullfile(release_root,'tmp');
if ~exist(out_dir,'dir'); mkdir(out_dir); end
writetable(results,fullfile(out_dir,'smoke_test_v5.csv'));
assert(all(strcmp(results.status,'PASS')), ...
    'HEBCPF:SmokeTest','At least one v5 smoke-test run failed.');
delete(gcp('nocreate'));
end

function row = empty_row()
row = struct('implementation','','case_name','','policy','','status','', ...
    'solutions',NaN,'expected_solutions',NaN,'max_residual',NaN, ...
    'wall_time_sec',NaN,'total_traces',NaN,'trace_to_90pct',NaN, ...
    'message','');
end

function value = ternary(condition,true_value,false_value)
if condition; value = true_value; else; value = false_value; end
end
