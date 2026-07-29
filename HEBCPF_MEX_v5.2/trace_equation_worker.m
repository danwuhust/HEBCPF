function [Z, Time, NY_inc] = trace_equation_worker(MS, M0, Tinv, Ma, Mp, Mq, Ybus, ...
    x, equationnumber, solutionnumber, param, bb, ba, km0, bus_n, I, degree, numofcons, svec)
%TRACE_EQUATION_WORKER  Run all-solution tracing for a single equation (called by parfeval)
%
%  Inputs:
%    MS, M0, Tinv, Ma, Mp, Mq, Ybus  - system matrices
%    x                - initial point of the current solution
%    equationnumber   - index of the equation being traced
%    solutionnumber   - index of the current solution
%    param            - parameter struct (contains fac, sstep, targetsteps, etc.)
%    bb, ba, km0      - equation vectors
%    bus_n, I, degree - system dimensions and holomorphic parameters
%    numofcons        - number of constraints
%    svec             - equation selection vector
%
%  Outputs:
%    Z       - all solutions found by this trace (numofcons x k)
%    Time    - trace time (seconds)
%    NY_inc  - increment in the number of columns of the Y matrix

    persistent workerThreadLimited
    if isempty(workerThreadLimited)
        maxNumCompThreads(1);
        workerThreadLimited = true;
    end
    warning off MATLAB:nearlySingularMatrix;
    Fail = -1;
    NYh = 0;           % used by script line 338: NYh=NYh+size(aal,1)

    tic;

    % call the tracing script (runs in this function's workspace, sharing all local variables)
    % hybrid_traceloops_4_with_mex.m requires the following variables in the workspace:
    %   x, equationnumber, solutionnumber, param, bb, ba, km0,
    %   MS, M0, Tinv, Ma, Mp, Mq, Ybus, bus_n, I, degree, numofcons, svec,
    %   Fail, NYh
    % after the script runs it produces: Z (solution matrix), Y (trace path)
    fac = param.fac;   % script line 25 param.fac=fac requires fac to exist
    script_dir = fileparts(mfilename('fullpath'));
    run(fullfile(script_dir, 'hybrid_traceloops_4_with_mex.m'));

    Time = toc;
    NY_inc = size(Y, 2);
end
