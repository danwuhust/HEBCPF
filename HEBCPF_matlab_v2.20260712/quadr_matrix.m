function [Ma,ra,npg,npd,nq,nv,I] = quadr_matrix(Mp1,Mq1,mpc,bus_n,gen_n,numofvar)
% global bus_n numofvar gen_n;
bus=mpc.bus;
gen=mpc.gen;
Ms=zeros(2*bus_n,2*bus_n);
Mpg=zeros(2*bus_n,2*bus_n);
Mpd=zeros(2*bus_n,2*bus_n);
Mq=zeros(2*bus_n,2*bus_n);
Mv=zeros(2*bus_n,2*bus_n);
npg=0;
npgt=numofvar;
npd=0;
nq=0;
nv=0;
ns=numofvar;
I.pv=[];
I.pq=[];
for i=1:1:bus_n
    % slack bus
    if bus(i,2)==3
       Ms(i,i)=1;
       Ms(bus_n+i,bus_n+i)=1;
%        rs=bus(i,8)^2;
       for j=1:gen_n
           if gen(j,1)==i
                rs=gen(j,6)^2;
           end
       end
       ns=i;
       npgt=npg+1;
    end
    % PV
    if bus(i,2)==2
       npg=npg+1;
       I.pv(npg,1)=i;
       Mpg(:,:,npg)=full(Mp1{i}); 
       if gen(npg,1)<ns
          rpg(npg,1)=(gen(npg,2)-bus(i,3))/mpc.baseMVA;
       else
          rpg(npg,1)=(gen(npg+1,2)-bus(i,3))/mpc.baseMVA; 
       end
       
       nv=nv+1;
       Mv(:,:,nv)=zeros(2*bus_n,2*bus_n);
       Mv(i,i,nv)=1;
       Mv(i+bus_n,i+bus_n,nv)=1;
       if nv<npgt
          rv(nv,1)=gen(nv,6)^2;
       else
          rv(nv,1)=gen(nv+1,6)^2; 
       end
    end
    % PQ
    if bus(i,2)==1
       npd=npd+1;
       I.pq(npd,1)=i;
       Mpd(:,:,npd)=full(Mp1{i});
       rpd(npd,1)=(-bus(i,3))/mpc.baseMVA;
        
       nq=nq+1;
       Mq(:,:,nq)=full(Mq1{i});
       rq(nq,1)=(-bus(i,4))/mpc.baseMVA;
    end
end
I.slack=ns;
%
Ms(bus_n+ns,:)=[];
Ms(:,bus_n+ns)=[];
Mpg(bus_n+ns,:,:)=[];
Mpg(:,bus_n+ns,:)=[];
Mpd(bus_n+ns,:,:)=[];
Mpd(:,bus_n+ns,:)=[];
Mq(bus_n+ns,:,:)=[];
Mq(:,bus_n+ns,:)=[];
Mv(bus_n+ns,:,:)=[];
Mv(:,bus_n+ns,:)=[];
%}
Ma=cell(2*bus_n-1,1);
Ma{1}=sparse(Ms);
ra(1,1)=rs;

for i=1:1:nv
    Ma{1+i}=sparse(Mv(:,:,i));
    ra(1+i,1)=rv(i,1);
end
for i=1:1:nq
    Ma{1+nv+i}=sparse(Mq(:,:,i));
    ra(1+nv+i,1)=rq(i,1);
end
for i=1:1:npg
    Ma{1+nv+nq+i}=sparse(Mpg(:,:,i));
    ra(1+nv+nq+i,1)=rpg(i,1);
end
for i=1:1:npd
    Ma{1+nv+nq+npg+i}=sparse(Mpd(:,:,i));
    ra(1+nv+nq+npg+i,1)=rpd(i,1);
end
%{
Ma(ns,ns,numofvar)=1;
Ma(bus_n+ns,bus_n+ns,numofvar)=1;
Ma(ns,bus_n+ns,numofvar)=0.5;
Ma(bus_n+ns,ns,numofvar)=0.5;
ra(end+1)=bus(ns,8)^2;
%}
end