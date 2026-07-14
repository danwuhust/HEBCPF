function [Vd,Vq,xpade,dbap]=update_V_Pade(num,den,degree,Td,holo_arc,dst,bus_n,I)
% dst is scalar, so evaluate each Padé polynomial by one matrix-vector product.
Vnum = num * (dst.^(0:degree.numerator)).';
Vden = den * (dst.^(0:degree.denominator)).';
V=Vnum./Vden;
Vd=real(V);
Vq=imag(V);

xpade=[Vd;Vq];
xpade(bus_n+I.slack)=[];
dbap=Td*dst*holo_arc;

end
