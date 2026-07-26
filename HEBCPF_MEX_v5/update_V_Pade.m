function [Vd,Vq,xpade,dbap]=update_V_Pade(num,den,degree,Td,holo_arc,dst,bus_n,I)
% [PROTO E] dst is a SCALAR: the stock code built two bus_n x (deg+1)
% matrices whose every column is the same scalar power, then did an
% elementwise multiply and a row sum. Equivalent (and allocation-free):
% one power vector and one BLAS mat-vec per polynomial.

Vnum = num * (dst.^(0:degree.numerator)).';
Vden = den * (dst.^(0:degree.denominator)).';
V=Vnum./Vden;
Vd=real(V);
Vq=imag(V);

xpade=[Vd;Vq];
xpade(bus_n+I.slack)=[];
dbap=Td*dst*holo_arc;

end
