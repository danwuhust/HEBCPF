function [Vx,Vy,Wx,Wy,Qpv,Kt]=holomorphic_cont_tri(Td,holo_arc,Ppm,Qqm,Lh,Uh,Start,I,degree,bus_n)

% global bus_n

Kt=Td*holo_arc;
% Flattened from local structs K/dK to plain variables so MATLAB Coder
% (codegen) does not see fields that are "undefined on some execution paths".
Ksl=Kt(1);
Kpvv=zeros(bus_n,1);
Kpvv(I.pv)=Kt(2:1+size(I.pv,1));
dKpvp=zeros(bus_n,1);
dKpvp(I.pv)=Kt(bus_n+1:bus_n+size(I.pv,1));
dKpqq=zeros(bus_n,1);
dKpqq(I.pq)=Kt(1+size(I.pv,1)+1:1+size(I.pv,1)+size(I.pq,1));
dKpqp=zeros(bus_n,1);
dKpqp(I.pq)=Kt(bus_n+size(I.pv,1)+1:end);

Rsl=zeros(1,degree.max);
% Rpv=zeros(size(I.pv,1),degree.max);
% Rpq=zeros(size(I.pq,1),degree.max);

Vx=zeros(bus_n,degree.max);
Vy=zeros(bus_n,degree.max);
Wx=zeros(bus_n,degree.max);
Wy=zeros(bus_n,degree.max);
Qpv=zeros(bus_n,degree.max);
Vx(:,1)=real(Start.V0);
Vy(:,1)=imag(Start.V0);
Wx(:,1)=real(Start.W0);
Wy(:,1)=imag(Start.W0);
Qpv(:,1)=Start.Q0;
nslk=[1:bus_n];
nslk(find(nslk==I.slack))=[];
% Use the absolute injection changes directly. Dividing by the base
% injections and multiplying back is algebraically redundant and produces
% NaN/Inf for zero-injection buses.
kpq=dKpqp-1i*dKpqq;
kpo=dKpvp;
Nr=5*size(I.pv,1)+4*size(I.pq,1)+1;
dvq2=zeros(Nr,1);
dvq1=zeros(Nr,1);
dvq=zeros(Nr,1);

% Triangular solve options for linsolve (LAPACK dtrsv path)
optsL.LT = true;  optsL.UT = false; optsL.UHESS = false;
optsL.SYM = false; optsL.POSDEF = false; optsL.RECT = false; optsL.TRANSA = false;
optsU.LT = false; optsU.UT = true;  optsU.UHESS = false;
optsU.SYM = false; optsU.POSDEF = false; optsU.RECT = false; optsU.TRANSA = false;

for ii=1:degree.act
    if ii==1
        %% slack bus
        Rsl(ii)=Ksl/(2*Vx(I.slack,1));

        %% PV bus
        rpva=kpo.*(Wx(:,ii)-1i*Wy(:,ii));
        rpvb=zeros(bus_n,1);
        rpvc=Kpvv;
        
        %% PQ bus
        rpqa=kpq.*(Wx(:,ii)-1i*Wy(:,ii));
        rpqb=zeros(bus_n,1);
        
        %% right hand side
        rpq1=real(rpqa(I.pq,:));
        rpq2=imag(rpqa(I.pq,:));
        rpq3=real(rpqb(I.pq,:));
        rpq4=imag(rpqb(I.pq,:));

        rpv1=real(rpva(I.pv,:));
        rpv2=imag(rpva(I.pv,:));
        rpv3=real(rpvb(I.pv,:));
        rpv4=imag(rpvb(I.pv,:));
        rpv5=real(rpvc(I.pv,:));
        
        R=[Rsl(ii);rpv1;rpq1;rpv2;rpq2;rpv3;rpq3;rpv4;rpq4;rpv5];
        
        %% next degree
        Rpm=R(Ppm);
        % solving triangle linear system (LAPACK dtrsv via linsolve)
        dvq1 = linsolve(Lh, Rpm, optsL);
        dvq2 = linsolve(Uh, dvq1, optsU);

        dvq=dvq2(Qqm);
        %
        Vx(:,ii+1)=dvq(1:bus_n);
        Vy(nslk,ii+1)=dvq(bus_n+1:2*bus_n-1);
        Wx(nslk,ii+1)=dvq(2*bus_n:3*bus_n-2);
        Wy(nslk,ii+1)=dvq(3*bus_n-1:4*bus_n-3);
        Qpv(I.pv,ii+1)=dvq(4*bus_n-2:end);

    else
        vcm=Vx(:,2:ii)+1i*Vy(:,2:ii);
        wcm=Wx(:,2:ii)+1i*Wy(:,2:ii);
        %% slack bus       
        Rsl(ii)=-Rsl(1:ii-1)*flip(Rsl(1:ii-1)')/(2*Vx(I.slack,1));
        
        %% PQ bus
        rpqa=kpq.*(Wx(:,ii)-1i*Wy(:,ii));
%         rpqb=-sum((Vx(:,2:ii)+1i*Vy(:,2:ii)).*flip((Wx(:,2:ii)+1i*Wy(:,2:ii)),2),2);
        rpqb=-sum(vcm.*flip(wcm,2),2);
        %% PV bus
%         rpva=kpo.*(Wx(:,ii)-1i*Wy(:,ii))-1i*(sum(Qpv(:,2:ii).*flip(Wx(:,2:ii)-1i*Wy(:,2:ii),2),2));
        rpva=kpo.*(Wx(:,ii)-1i*Wy(:,ii))-1i*(sum(Qpv(:,2:ii).*flip(conj(wcm),2),2));
%         rpvb=-sum((Vx(:,2:ii)+1i*Vy(:,2:ii)).*flip((Wx(:,2:ii)+1i*Wy(:,2:ii)),2),2);
        rpvb=-sum(vcm.*flip(wcm,2),2);
%         rpvc=-sum((Vx(:,2:ii)+1i*Vy(:,2:ii)).*flip((Vx(:,2:ii)-1i*Vy(:,2:ii)),2),2);
        rpvc=-sum(vcm.*flip(conj(vcm),2),2);
        
        %% right hand side
        rpq1=real(rpqa(I.pq,:));
        rpq2=imag(rpqa(I.pq,:));
        rpq3=real(rpqb(I.pq,:));
        rpq4=imag(rpqb(I.pq,:));

        rpv1=real(rpva(I.pv,:));
        rpv2=imag(rpva(I.pv,:));
        rpv3=real(rpvb(I.pv,:));
        rpv4=imag(rpvb(I.pv,:));
        rpv5=real(rpvc(I.pv,:));
        
        R=[Rsl(ii);rpv1;rpq1;rpv2;rpq2;rpv3;rpq3;rpv4;rpq4;rpv5];
        
        %% next degree
        Rpm=R(Ppm);
        % solving triangle linear system (LAPACK dtrsv via linsolve)
        dvq1 = linsolve(Lh, Rpm, optsL);
        dvq2 = linsolve(Uh, dvq1, optsU);

        dvq=dvq2(Qqm);
        %

        Vx(:,ii+1)=real(dvq(1:bus_n));
        Vy(nslk,ii+1)=real(dvq(bus_n+1:2*bus_n-1));
        Wx(nslk,ii+1)=real(dvq(2*bus_n:3*bus_n-2));
        Wy(nslk,ii+1)=real(dvq(3*bus_n-1:4*bus_n-3));
        Qpv(I.pv,ii+1)=real(dvq(4*bus_n-2:end));
    end
end

end
