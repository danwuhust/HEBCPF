function [solu,solu_0] = Solu(result,ns,bus_n)
% global bus_n
vm=result.bus(:,8);
angle=result.bus(:,9)*pi/180;
vd=vm.*cos(angle);
vq=vm.*sin(angle);

solu_0=[vd;vq];
solu=solu_0;
solu(bus_n+ns)=[];
end