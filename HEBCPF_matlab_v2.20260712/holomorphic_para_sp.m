function [H,Start]=holomorphic_para_sp(Mp,Mq,solu0,Ybus,I,bus_n)
% global bus_n
P0=zeros(bus_n,1);
Q0=zeros(bus_n,1);
V0=solu0(1:bus_n)+1i*solu0(bus_n+1:end);
W0=1./V0;
for i=1:bus_n
    P0(i,1)=solu0'*Mp{i}*solu0;
    Q0(i,1)=solu0'*Mq{i}*solu0;
end
G=real(Ybus);
B=imag(Ybus);
Ap=diag(P0);
Aq=diag(Q0);
Avx=sparse(diag(real(V0)));
Avy=sparse(diag(imag(V0)));
Awx=sparse(diag(real(W0)));
Awy=sparse(diag(imag(W0)));
spz=sparse(zeros(bus_n,bus_n));

Mpq=[G -B -Ap Aq spz; 
    B G Aq Ap spz; 
    Awx -Awy Avx -Avy spz; 
    Awy Awx Avy Avx spz];
Mpv=[G -B -Ap Aq Awy; 
    B G Aq Ap Awx; 
    Awx -Awy Avx -Avy spz; 
    Awy Awx Avy Avx spz;
    2*Avx 2*Avy spz spz spz];
dl=[bus_n+I.slack 2*bus_n+I.slack 3*bus_n+I.slack 4*bus_n+I.pq' 4*bus_n+I.slack];
Mpq(:,dl)=[];
Mpv(:,dl)=[];

Hpq1=Mpq(I.pq,:);
Hpq2=Mpq(bus_n+I.pq,:);
Hpq3=Mpq(2*bus_n+I.pq,:);
Hpq4=Mpq(3*bus_n+I.pq,:);

Hpv1=Mpv(I.pv,:);
Hpv2=Mpv(bus_n+I.pv,:);
Hpv3=Mpv(2*bus_n+I.pv,:);
Hpv4=Mpv(3*bus_n+I.pv,:);
Hpv5=Mpv(4*bus_n+I.pv,:);

Hsl=zeros(1,size(Mpv,2));
Hsl(I.slack)=1;

H=[Hsl;Hpv1;Hpq1;Hpv2;Hpq2;Hpv3;Hpq3;Hpv4;Hpq4;Hpv5];

Start.V0=V0;
Start.W0=W0;
Start.P0=P0;
Start.Q0=Q0;

end