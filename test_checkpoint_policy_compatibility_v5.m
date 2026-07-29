function results = test_checkpoint_policy_compatibility_v5()
%TEST_CHECKPOINT_POLICY_COMPATIBILITY_V5 Cross-policy resume regression.
% Creates a partial case5loop checkpoint under each implementation and resumes
% it with a different canonical scheduler. The forced legacy time threshold is
% test-only; production checkpoint cadence is unchanged.

release_root = fileparts(mfilename('fullpath'));
implementations = {'HEBCPF_MEX_v5.2','HEBCPF_matlab_v5.2'};
transitions = {'scan','bandit'; 'bandit','novelty'; 'novelty','scan'};
rows = repmat(empty_row(),0,1);

p = gcp('nocreate');
if isempty(p); p = parpool('local'); end
for ii = 1:numel(implementations)
    solver_dir = fullfile(release_root,implementations{ii});
    for ti = 1:size(transitions,1)
        source_policy = transitions{ti,1};
        resume_policy = transitions{ti,2};
        checkpoint_file = create_checkpoint(solver_dir,source_policy);
        checkpoint = load(checkpoint_file,'trace_policy','count','numberofsolutions');
        reference = resume_checkpoint(solver_dir,checkpoint_file,source_policy);
        resumed = resume_checkpoint(solver_dir,checkpoint_file,resume_policy);
        row = empty_row();
        row.implementation = implementations{ii};
        row.source_policy = source_policy;
        row.resume_policy = resume_policy;
        row.checkpoint_policy = checkpoint.trace_policy;
        row.checkpoint_traces = double(checkpoint.count);
        row.checkpoint_solutions = double(checkpoint.numberofsolutions);
        row.final_solutions = size(resumed,2);
        if isequal(size(reference),size(resumed))
            row.max_abs_difference = max(abs(reference(:)-resumed(:)));
        else
            row.max_abs_difference = Inf;
        end
        row.status = ternary(isfinite(row.max_abs_difference) && ...
            row.max_abs_difference <= 4e-7,'PASS','FAIL');
        rows(end+1) = row; %#ok<AGROW>
        fprintf('%s %s -> %s: %s (%d checkpoint solutions, max diff %.3g)\n', ...
            implementations{ii},source_policy,resume_policy,row.status, ...
            row.checkpoint_solutions,row.max_abs_difference);
    end
    checkpoint_file = fullfile(solver_dir,'temp_result.mat');
    test_checkpoint_file = fullfile(solver_dir,'temp_result_policy_test.mat');
    if exist(checkpoint_file,'file'); delete(checkpoint_file); end
    if exist(test_checkpoint_file,'file'); delete(test_checkpoint_file); end
end
results = struct2table(rows,'AsArray',true);
out_dir = fullfile(release_root,'tmp');
if ~exist(out_dir,'dir'); mkdir(out_dir); end
writetable(results,fullfile(out_dir,'checkpoint_policy_compatibility_v5.csv'));
assert(all(strcmp(results.status,'PASS')), ...
    'HEBCPF:CheckpointCompatibility','At least one cross-policy resume failed.');
end

function checkpoint_file = create_checkpoint(solver_dir,policy)
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir)); %#ok<NASGU>
cd(solver_dir);
addpath(solver_dir);
global HEBCPOLICY %#ok<GVMIS>
HEBCPOLICY = policy;
source = fileread('main.m');
active = regexp(source,'^[ \t]*mpc=case\w+;','match','once','lineanchors');
source = strrep(source,active,'mpc=case5loop;');
evalc('eval(source);');
setenv('HEBCPF_CHECKPOINT_TRACE_INTERVAL','3');
setenv('HEBCPF_STOP_AFTER_CHECKPOINT','true');
cleanup_env = onCleanup(@() clear_checkpoint_environment()); %#ok<NASGU>
evalc('run(''runVBook_hybrid_parfeval.m'');');
checkpoint_file = fullfile(pwd,'temp_result_policy_test.mat');
copyfile(fullfile(pwd,'temp_result.mat'),checkpoint_file,'f');
end

function final_solutions = resume_checkpoint(solver_dir,checkpoint_file,policy)
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir)); %#ok<NASGU>
cd(solver_dir);
addpath(solver_dir);
load(checkpoint_file);
global HEBCPOLICY %#ok<GVMIS>
HEBCPOLICY = policy;
evalc('run(''runVBook_hybrid_parfeval.m'');');
final_solutions = Zsave;
end

function clear_checkpoint_environment()
setenv('HEBCPF_CHECKPOINT_TRACE_INTERVAL','');
setenv('HEBCPF_STOP_AFTER_CHECKPOINT','');
end

function row = empty_row()
row = struct('implementation','','source_policy','','resume_policy','', ...
    'checkpoint_policy','','checkpoint_traces',NaN, ...
    'checkpoint_solutions',NaN,'final_solutions',NaN, ...
    'max_abs_difference',NaN,'status','');
end

function value = ternary(condition,yes,no)
if condition; value = yes; else; value = no; end
end
