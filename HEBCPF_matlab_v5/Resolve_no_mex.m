function [x,flag] = Resolve_no_mex(MS,bb,M0,km0,x,n,ind,param)
% Newton corrector. Constants come from solver_params via the optional `param`;
% if omitted (e.g. the one-time initial resolve before param exists) the defaults
% below are used and are identical to solver_params' resolve_* values.
% Pure-MATLAB v4 deliberately uses MATLAB's built-in linear algebra only.
if nargin>=8 && isfield(param,'resolve_maxiter')
    maxiter=param.resolve_maxiter; err=param.resolve_err;
    err1=param.resolve_err1;       fac=param.resolve_fac;
else
    maxiter=10; err=2e-10; err1=1e-7; fac=0.4e1;
end
flag=0;
for kk=1:1:maxiter
    d0=x'*M0*x;
    b=2*M0*x;
    Jac = MakeJacobianD(x,MS,n);
    g   = Jac*x/2 + d0*km0 - bb;
    J   = [Jac fac*km0; b'/fac^2 -1/fac];
    F   = [g;0];
    dy  = -J\F;
    dx=dy(1:end-1);
    x=x+dx;
    if ind==1
        fprintf('No. Iter: %d, Step length: %g, Norm of gap: %g \n',kk,norm(dx),norm(g));
    end
    if (norm(g)<=err && norm(dx)<=err1)
        flag=1;
        if ind==1
            fprintf('Initial guess resolved! \n');
        end
       break;
    end
end
if kk == maxiter
    flag=-1;
    if ind==1
        fprintf('Exceed maximun iteration limit! \n');
    end
end
end
