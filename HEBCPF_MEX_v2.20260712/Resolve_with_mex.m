function [x,flag] = Resolve_with_mex(MS,bb,M0,km0,x,n,ind)
% global gen_n1 bus_n1 branch_n1;
% [KLU] The bordered Newton solve uses the KLU symbolic-once refactor (klurf) when the
% MEX is present (dense fallback otherwise) - same auto default as branch_trace; set
% global USEKLU to override (0 dense, 1 klu full-factor, 2 klurf refactor).
persistent klurf_ok
if isempty(klurf_ok); klurf_ok = (exist('klurf','file')==3); end
global USEKLU; %#ok<GVMIS>
maxiter=10;
err=2e-10;
err1=1e-7;
fac=0.4e1;
flag=0;
if isempty(USEKLU); klumode = 2*klurf_ok; else; klumode = USEKLU; end
for kk=1:1:maxiter
    d0=x'*M0*x;
    b=2*M0*x;
    if klumode>0
        % bordered system [A, fac*km0; b'/fac^2, -1/fac], A = sparse power-flow core.
        Jc = MakeJacobianD(x,MS,n,'sparse');
        g  = Jc*x/2 + d0*km0 - bb;
        if klumode==2
            GH = klurf(Jc, [g, fac*km0]);              % analyze-once then refactor
        else
            LUk = klu(Jc); GH = klu(LUk,'\',[g, fac*km0]);
        end
        gk = GH(:,1); Hk = GH(:,2);
        Ck = b'/fac^2;
        y2 = (-1/fac - Ck*Hk) \ (-(Ck*gk));            % 1x1 Schur (f2 = 0)
        dy = -[gk - Hk*y2; y2];
    else
        Jac = MakeJacobianD(x,MS,n);
        g   = Jac*x/2 + d0*km0 - bb;
        J   = [Jac fac*km0; b'/fac^2 -1/fac];
        F   = [g;0];
        dy  = -J\F;
    end
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
