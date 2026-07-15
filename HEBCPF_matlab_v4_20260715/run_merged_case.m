function [result,solutions] = run_merged_case(case_name,mode,close_pool_on_finish)
%RUN_MERGED_CASE Execute one case with serial, parfor, or parfeval tracing.
%
% close_pool_on_finish defaults to true for non-desktop MATLAB sessions,
% including "matlab -batch ...".  This avoids a Windows/R2022a shutdown
% crash where MATLAB exits while worker payloads are still being destroyed.

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
    result=struct('case_name',getenv('MERGED_CASE_NAME'), ...
        'mode',getenv('MERGED_SOLVER_MODE'),'status','PASS', ...
        'buses',bus_n,'solutions',size(Zsave,2), ...
        'max_residual',max_residual,'wall_time',wall_time, ...
        'error_id','','message','');
catch exception
    solutions=[];
    result=struct('case_name',getenv('MERGED_CASE_NAME'), ...
        'mode',getenv('MERGED_SOLVER_MODE'),'status','FAIL', ...
        'buses',NaN,'solutions',NaN,'max_residual',NaN, ...
        'wall_time',NaN,'error_id',exception.identifier, ...
        'message',regexprep(exception.message,'[\r\n]+',' '));
end
clear futures MS_c M0_c Tinv_c Ma_c Mp_c Mq_c Ybus_c
if strcmp(getenv('MERGED_SOLVER_MODE'),'parfeval') || strcmp(getenv('MERGED_SOLVER_MODE'),'parfor')
    cleanup_parallel_artifacts(strcmp(getenv('MERGED_CLOSE_POOL_ON_FINISH'),'true'));
end
end
