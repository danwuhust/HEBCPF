function cleanup_parallel_artifacts(close_pool)
%CLEANUP_PARALLEL_ARTIFACTS Release parallel payloads before MATLAB exits.
%
% MATLAB R2022a on Windows can crash during -batch shutdown if a local
% parpool is still holding large parallel.pool.Constant payloads or recently
% completed FevalFuture objects.  This helper clears solver-specific worker
% caches and can optionally delete the pool.

if nargin<1
    close_pool=false;
end

if license('test','Distrib_Computing_Toolbox')
    pool=gcp('nocreate');
    if ~isempty(pool)
        try
            pctRunOnAll clear MakeJacobianD quad_form_vals trace_equation_worker trace_equation_worker_const
        catch
        end
        if close_pool
            try
                delete(pool);
            catch
            end
        end
    end
end

clear MakeJacobianD quad_form_vals trace_equation_worker trace_equation_worker_const
end
