function [fac] = Tune_Fac(MS,M0,km0,x,n)
% fac=0;
% d0=x'*M0*x;
% Jac=[];
% J=[];
% F=[];
re0=1e16;
for fack=0:1:31
    fact=fack*2+1;
%     for i=1:1:size(bb)
%         Jac(i,:)=2*x'*full(M{i});
%     end    
%     b0=2*full(M0)*x;
%     J=[Jac fact*km0;b0'/fact^2 -1/fact];
%     eJ=abs(eig(J'*J));
%     re1=sqrt(max(eJ)/min(eJ));
%     [SJ] = MakeJacobian_mex(x,MS,n,1);
    [SJ] = MakeJacobian(x,MS,n,1);
    Jac=2*sparse(SJ(:,2),SJ(:,1),SJ(:,3),n,n)';    
    b0=2*M0*x;
    J=[Jac fact*km0;b0'/fact^2 -1/fact];
    eJ=abs(eig(J'*J));
    re1=sqrt(max(eJ)/min(eJ));
    fprintf('Scaling factor: %g; Condition number of Jacobian matrix: %g\n',fact,re1);
    if re1<re0        
       re0=re1;
       fac=fact;
    end
end

end