function [H,Start]=holomorphic_para_sp(Mp,Mq,solu0,Ybus,I,bus_n)
% Vectorized quadratic forms and cached case invariants.

persistent VP VQ Gc Bc spzc dlc idxc nc nnzYc sumYc
key_n = bus_n; key_nnz = nnz(Ybus); key_sum = full(sum(abs(nonzeros(Ybus))));
if isempty(nc) || key_n~=nc || key_nnz~=nnzYc || key_sum~=sumYc
    N = 2*bus_n;
    rows = cell(2*bus_n,1); cols = cell(2*bus_n,1); vals = cell(2*bus_n,1);
    for i = 1:bus_n
        [r,c,v] = find(Mp{i});
        rows{i} = r + (i-1)*N; cols{i} = c; vals{i} = v;
        [r,c,v] = find(Mq{i});
        rows{bus_n+i} = r + (bus_n+i-1)*N; cols{bus_n+i} = c; vals{bus_n+i} = v;
    end
    VPQ = sparse(cell2mat(rows), cell2mat(cols), cell2mat(vals), 2*bus_n*N, N);
    VP = VPQ(1:bus_n*N,:); VQ = VPQ(bus_n*N+1:end,:);
    Gc = real(Ybus); Bc = imag(Ybus); spzc = sparse(bus_n,bus_n);
    dlc = [bus_n+I.slack 2*bus_n+I.slack 3*bus_n+I.slack 4*bus_n+I.pq' 4*bus_n+I.slack];
    idxc = (1:bus_n)'; nc = key_n; nnzYc = key_nnz; sumYc = key_sum;
end

V0 = solu0(1:bus_n) + 1i*solu0(bus_n+1:end);
W0 = 1./V0;
N = 2*bus_n;
P0 = reshape(VP*solu0, N, bus_n).' * solu0;
Q0 = reshape(VQ*solu0, N, bus_n).' * solu0;
G=Gc; B=Bc; spz=spzc;
Ap=sparse(idxc,idxc,P0,bus_n,bus_n);
Aq=sparse(idxc,idxc,Q0,bus_n,bus_n);
Avx=sparse(idxc,idxc,real(V0),bus_n,bus_n);
Avy=sparse(idxc,idxc,imag(V0),bus_n,bus_n);
Awx=sparse(idxc,idxc,real(W0),bus_n,bus_n);
Awy=sparse(idxc,idxc,imag(W0),bus_n,bus_n);

Mpq=[G -B -Ap Aq spz; 
    B G Aq Ap spz; 
    Awx -Awy Avx -Avy spz; 
    Awy Awx Avy Avx spz];
Mpv=[G -B -Ap Aq Awy; 
    B G Aq Ap Awx; 
    Awx -Awy Avx -Avy spz; 
    Awy Awx Avy Avx spz;
    2*Avx 2*Avy spz spz spz];
Mpq(:,dlc)=[];
Mpv(:,dlc)=[];

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
