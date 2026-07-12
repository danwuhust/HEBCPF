function [SJ] = MakeJacobian(b,S,n,dd)
% generate Jacobian matrix
N=size(S,1);
SJ=zeros(N,3);
tol=1e-14;
% tol2=1e-15;
if dd==0
    for i=1:N
        a1=fix(S(i,2)/n);
        a2=mod(S(i,2),n);
        if a2==0
            a2=n;
            a1=a1-1;
        end
        SJ(i,1)=a1+1;
        SJ(i,2)=a2;
        SJ(i,3)=SJ(i,3)+S(i,3)*b(S(i,1));
        if abs(SJ(i,3))<tol
            SJ(i,3)=0;
        end
    end
else
    for i=1:N
        a1=fix(S(i,2)/n);
        a2=mod(S(i,2),n);
        if a2==0
            a2=n;
            a1=a1-1;
        end
        SJ(i,1)=a1+1;
        SJ(i,2)=a2;
        SJ(i,3)=SJ(i,3)+S(i,3)*b(S(i,1));
%         if abs(SJ(i,3))<tol2
%             SJ(i,3)=0;
%         end
    end
end

end
