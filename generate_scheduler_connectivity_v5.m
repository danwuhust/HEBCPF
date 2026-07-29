function generate_scheduler_connectivity_v5()
%GENERATE_SCHEDULER_CONNECTIVITY_V5 Case14mod diagrams for all queue policies.
release_root = fileparts(mfilename('fullpath'));
solver_dir = fullfile(release_root,'HEBCPF_MEX_v5.2');
policies = {'scan','bandit','novelty'};
cases = {'case14mod','case39'};
p = gcp('nocreate');
if isempty(p); p = parpool('local'); end %#ok<NASGU>
for ci = 1:numel(cases)
    for pi = 1:numel(policies)
        generate_one(solver_dir,cases{ci},policies{pi});
    end
end
end

function generate_one(solver_dir,case_name,policy)
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir)); %#ok<NASGU>
cd(solver_dir); addpath(solver_dir);
global HEBCPOLICY %#ok<GVMIS>
HEBCPOLICY = policy;
setenv('HEBCPF_CONNECTIVITY_CASE',case_name);
setenv('HEBCPF_CONNECTIVITY_ROOT',fileparts(solver_dir));
source = fileread('main.m');
active = regexp(source,'^[ \t]*mpc=case\w+;','match','once','lineanchors');
source = strrep(source,active,'mpc=feval(getenv(''HEBCPF_CONNECTIVITY_CASE''));');
evalc('eval(source);');
evalc('run(''runVBook_hybrid_parfeval.m'');');
release_root = getenv('HEBCPF_CONNECTIVITY_ROOT');
case_name = getenv('HEBCPF_CONNECTIVITY_CASE');
out_png = fullfile(release_root,sprintf('connectivity_%s_%s.png',case_name,trace_policy));
solution_connectivity(VBook,out_png,'Dim',2,'Restarts',8,'Start',solu,'Zsave',Zsave);
end
