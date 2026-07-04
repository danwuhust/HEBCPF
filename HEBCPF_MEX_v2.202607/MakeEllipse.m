function [MS,M,bb,M0,km0,bm0,succ,Mss,T] = MakeEllipse(Ma, b,maxtrial,bus_n,gen_n,numofvar,numofcons)
% global bus_n gen_n numofvar numofcons;
Min=0.001;
% Mm=zeros(numofvandc,numofvandc,numofvandc);
borig=b;
bb=zeros(numofcons,1);
M0=zeros(numofvar,numofvar);
b0=0;
M1=zeros(numofvar,numofvar);
b1=0;
% for i=1:1:gen_n1
% %     M0 = M0 + Ma(:,:,i) - Ma(:,:,gen_n1+i) + Ma(:,:,2*gen_n1+i) - Ma(:,:,3*gen_n1+i);
%     M0 = M0 + full(Ma{i}) - full(Ma{gen_n1+i}) + full(Ma{2*gen_n1+i}) - full(Ma{3*gen_n1+i});
%     b0 = b0 + b(i) - b(gen_n1+i) + b(2*gen_n1+i) - b(3*gen_n1+i);
% end

Ta=eye(numofcons);
T=zeros(numofcons,numofcons);
LC0=zeros(1,numofcons);
for i=1:1:bus_n-gen_n
    M0 = M0 + Ma{gen_n+i};
    b0 = b0 + b(gen_n+i);
    LC0(gen_n+i)=1;
end

LC1=zeros(1,numofcons);
for i=1:1:gen_n
    M1 = M1 + Ma{i};
    b1 = b1 + b(i);
    LC1(i)=1;
end

k=0;
fk=1;
D = zeros(numofvar,1);
while (b0<0 || D(1,1)<Min) && k<maxtrial
    k=k+1;
    M0 = M0 + fk*k*M1;
    b0 = b0 + fk*k*b1;
    D=eig(M0);
    LC0=LC0+fk*k*LC1;
    fprintf('The %d-th trial to generate the original base ellipse. Scaling factor: %g, Minimum eigenvalue: %g \n',k,fk*k,D(1,1));
end

if k<maxtrial
succ=1;
% combine all the constraints into the base ellipse
%
dLC=ones(1,numofcons);
dM=zeros(numofvar,numofvar);
db=0;
for i=1:1:numofcons
    %{
    if abs(b(i))>0
        dM=dM+Ma{i}/b(i);
        db=db+1;
    else
        dM=dM+Ma{i};
        db=db+b(i);
    end
    %}
    %
    sm=max(max(abs(Ma{i})));
    Ma{i}=Ma{i}/sm;
    b(i)=b(i)/sm;
    dLC(i)=dLC(i)/sm;
    Ta(i,i)=Ta(i,i)/sm;
    %}
    dM=dM+Ma{i};
    db=db+b(i);
end
%}
tb=0;
td=zeros(numofvar,1);
tk=0;
ftk=2;
while (tb<0 || td(1)<Min) 
    tk=tk+1;
    tM = dM + ftk*tk*M0;
    tb = db + ftk*tk*b0;
    td=eig(tM);
    tLC = dLC + ftk*tk*LC0;
    fprintf('The %d-th trial to generate the actual base ellipse. Scaling factor: %g \n',tk,ftk*tk);
end
M0=tM;
b0=tb;
LC0 = tLC;
%}
fprintf('The base ellipse generated! \n');

% DD=zeros(numofvar,numofvar);
% [U,D]=eig(M0);
% for i=1:1:numofvar
%     DD(i,i)=sqrt(1/D(i,i));
% end
% T=U'\DD;

d = zeros(numofvar,numofvar);
ks = zeros(numofcons,1);
M=cell(1);
fa=5;
tol2=1e-4;
bm0=zeros(size(bb,1),1);
for i=1:1:numofcons
    %
    if abs(b(i))>tol2
       Ta(i,i)=Ta(i,i)/abs(b(i));
       Ma{i}=Ma{i}/abs(b(i));
       b(i)=b(i)/abs(b(i));
    end
    %}
    %{
    sm=max(max(abs(Ma{i})));
    Ma{i}=Ma{i}/sm;
    b(i)=b(i)/sm;
    %}
    while (bm0(i)<0 || d(1,i)<Min) 
        ks(i)=ks(i)+1;
        Mm=Ma{i}+fa*ks(i)*M0;
%         Mm(:,:,i)=T'*(Ma(:,:,i)+2*k(i)*M0)*T;
%         Mm(:,:,i)=0.5*(Mm(:,:,i)'+Mm(:,:,i));
        bm0(i,1)=b(i)+fa*ks(i)*b0;
        tm=Ta(i,:)+fa*ks(i)*LC0;
        d(:,i)=eig(Mm);  
    end
    M{i}=Ma{i}/bm0(i);
    bb(i)=1;
    T(i,:)=tm/bm0(i);
    fprintf('The %d-th ellise is generated! No. of trial %d. Scale factor %d. \n',i,ks(i),fa*ks(i));
end

M0=sparse(M0);
km0=fa*ks./bm0;
% Msp=cell(numofcons,1);
% for i=1:1:numofcons
%     [Msp{i}(:,1),Msp{i}(:,2),Msp{i}(:,3)]=find(M{i});
% end
Mss=cell2mat(M);
[MS(:,1),MS(:,2),MS(:,3)]=find(Mss);
else
    M=0;bb=0;M0=0;km0=0;
    bm0=0;succ=0;MS=0;Mss=0;
end

end