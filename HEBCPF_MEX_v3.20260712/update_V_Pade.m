function [Vd,Vq,xpade,dbap]=update_V_Pade(num,den,degree,Td,holo_arc,dst,bus_n,I)

numd=zeros(bus_n,degree.numerator+1);
for i=1:degree.numerator+1
    numd(:,i)=dst^(i-1);
end

dend=zeros(bus_n,degree.denominator+1);
for i=1:degree.denominator+1
    dend(:,i)=dst^(i-1);
end

Vnum=sum(num.*numd,2);
Vden=sum(den.*dend,2);
V=Vnum./Vden;
Vd=real(V);
Vq=imag(V);

xpade=[Vd;Vq];
xpade(bus_n+I.slack)=[];
dbap=Td*dst*holo_arc;

end