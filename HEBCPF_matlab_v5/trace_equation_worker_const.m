function [Z, Time, NY_inc] = trace_equation_worker_const(MS_c, M0_c, Tinv_c, Ma_c, Mp_c, Mq_c, Ybus_c, ...
    x, equationnumber, solutionnumber, param, bb, ba, km0, bus_n, I, degree, numofcons, svec)
%TRACE_EQUATION_WORKER_CONST  parfeval-specific wrapper
%  receives parallel.pool.Constant objects, unpacks them and calls trace_equation_worker

    MS   = MS_c.Value;
    M0   = M0_c.Value;
    Tinv = Tinv_c.Value;
    Ma   = Ma_c.Value;
    Mp   = Mp_c.Value;
    Mq   = Mq_c.Value;
    Ybus = Ybus_c.Value;

    [Z, Time, NY_inc] = trace_equation_worker(MS, M0, Tinv, Ma, Mp, Mq, Ybus, ...
        x, equationnumber, solutionnumber, param, bb, ba, km0, bus_n, I, degree, numofcons, svec);
end
