function [result,solutions] = run_merged_case(case_name,mode,close_pool_on_finish)
%RUN_MERGED_CASE Execute one case with serial, parfor, or parfeval tracing.
%
% close_pool_on_finish defaults to true for non-desktop MATLAB sessions,
% including "matlab -batch ...".  This avoids a Windows/R2022a shutdown
% crash where MATLAB exits while worker payloads are still being destroyed.

global HEBCPOLICY %#ok<GVMIS>
[HEBCPOLICY,~] = trace_policy_config(HEBCPOLICY);

setenv('MERGED_CASE_NAME',case_name);
setenv('MERGED_SOLVER_MODE',mode);
if nargin<3
    close_pool_on_finish=~usejava('desktop');
end
if close_pool_on_finish
    setenv('MERGED_CLOSE_POOL_ON_FINISH','true');
else
    setenv('MERGED_CLOSE_POOL_ON_FINISH','false');
end
try
    source=fileread('main.m');
    active=regexp(source,'^[ \t]*mpc=case\w+;','match','once','lineanchors');
    if isempty(active)
        error('MergedSolver:NoActiveCase','No active mpc=case... line in main.m.');
    end
    source=strrep(source,active,'mpc=feval(getenv(''MERGED_CASE_NAME''));');
    evalc('eval(source);');

    solve_clock=tic;
    switch getenv('MERGED_SOLVER_MODE')
        case 'serial'
            evalc('run(''runVBook_hybrid_2023.m'');');
        case 'parfor'
            evalc('run(''runVBook_hybrid_parallel.m'');');
        case 'parfeval'
            evalc('run(''runVBook_hybrid_parfeval.m'');');
        otherwise
            error('MergedSolver:UnknownMode','Unknown mode: %s', ...
                getenv('MERGED_SOLVER_MODE'));
    end
    wall_time=toc(solve_clock);

    max_residual=0;
    for solution_index=1:size(Zsave,2)
        z=Zsave(:,solution_index);
        for equation_index=1:numel(Ma)
            max_residual=max(max_residual, ...
                abs(z'*Ma{equation_index}*z-ba(equation_index)));
        end
    end
    solutions=Zsave;
    trace_summary = local_trace_summary();
    result=struct('case_name',getenv('MERGED_CASE_NAME'), ...
        'mode',getenv('MERGED_SOLVER_MODE'),'status','PASS', ...
        'buses',bus_n,'solutions',size(Zsave,2), ...
        'max_residual',max_residual,'wall_time',wall_time, ...
        'scheduler',HEBCPOLICY,'trace_metrics',trace_summary, ...
        'error_id','','message','');
catch exception
    solutions=[];
    result=struct('case_name',getenv('MERGED_CASE_NAME'), ...
        'mode',getenv('MERGED_SOLVER_MODE'),'status','FAIL', ...
        'buses',NaN,'solutions',NaN,'max_residual',NaN, ...
        'wall_time',NaN,'scheduler',HEBCPOLICY,'trace_metrics',struct(), ...
        'error_id',exception.identifier, ...
        'message',regexprep(exception.message,'[\r\n]+',' '));
end
clear futures MS_c M0_c Tinv_c Ma_c Mp_c Mq_c Ybus_c
if strcmp(getenv('MERGED_SOLVER_MODE'),'parfeval') || strcmp(getenv('MERGED_SOLVER_MODE'),'parfor')
    cleanup_parallel_artifacts(strcmp(getenv('MERGED_CLOSE_POOL_ON_FINISH'),'true'));
end
end

function summary = local_trace_summary()
vars = evalin('caller','whos');
names = {vars.name};
summary = struct('total_traces',NaN,'trace_to_50pct',NaN, ...
    'trace_to_90pct',NaN,'trace_to_99pct',NaN,'wall_to_90pct',NaN, ...
    'discovery_traces',NaN,'mean_worker_trace_sec',NaN, ...
    'p90_worker_trace_sec',NaN,'selection_wall_sec',NaN, ...
    'cursor_checks',NaN,'novelty_draws',NaN, ...
    'solution_history',[],'new_solutions_history',[], ...
    'worker_seconds_history',[],'completion_elapsed_history',[], ...
    'equation_history',[]);
if ~ismember('trace_solution_history',names); return; end
sol = double(evalin('caller','trace_solution_history(:)'));
worker = double(evalin('caller','trace_worker_seconds(:)'));
elapsed = double(evalin('caller','trace_completion_elapsed(:)'));
newsol = double(evalin('caller','trace_new_solutions(:)'));
if isempty(sol); return; end
summary.solution_history = sol;
summary.new_solutions_history = newsol;
summary.worker_seconds_history = worker;
summary.completion_elapsed_history = elapsed;
if ismember('trace_equation_history',names)
    summary.equation_history = double(evalin('caller','trace_equation_history(:)'));
end
final_n = sol(end);
offset = 0;
if ismember('trace_history_offset',names); offset = evalin('caller','trace_history_offset'); end
summary.total_traces = offset + numel(sol);
i50 = find(sol >= ceil(0.50*final_n),1,'first');
i90 = find(sol >= ceil(0.90*final_n),1,'first');
i99 = find(sol >= ceil(0.99*final_n),1,'first');
if ~isempty(i50); summary.trace_to_50pct = offset + i50; end
if ~isempty(i90); summary.trace_to_90pct = offset + i90; summary.wall_to_90pct = elapsed(i90); end
if ~isempty(i99); summary.trace_to_99pct = offset + i99; end
summary.discovery_traces = nnz(newsol > 0);
summary.mean_worker_trace_sec = mean(worker);
worker_sorted = sort(worker);
if ~isempty(worker_sorted)
    summary.p90_worker_trace_sec = worker_sorted(max(1,ceil(0.90*numel(worker_sorted))));
end
if ismember('scheduler_select_seconds',names)
    summary.selection_wall_sec = evalin('caller','scheduler_select_seconds');
end
if ismember('scheduler_cursor_checks',names)
    summary.cursor_checks = evalin('caller','scheduler_cursor_checks');
end
if ismember('scheduler_novelty_draws',names)
    summary.novelty_draws = evalin('caller','scheduler_novelty_draws');
end
end
