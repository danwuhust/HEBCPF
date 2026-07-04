function [xholo,aal1,aal_limit,aal_ind,holo_arc1,arc_ratio,num,den,kdst2,flagh]=...
    holomorphic_hybrid_5_no_mex(Mp,Mq,solu0,Ybus,I,Ma,bat,Tinv,holo_arc,holo_arcp,degree,...
    param,aal,aal_limit,aal_ind,sgn,bus_n,equationnumber,Display)

arc_ratio=sgn;
flagh=0;

[H,Start]=holomorphic_para_sp(Mp,Mq,solu0,Ybus,I,bus_n);

[Lh,Uh,Ph,Qh]=lu(H);
[Ppm,~,~]=find(Ph');
[Qqm,~,~]=find(Qh');
Lf=full(Lh);
Uf=full(Uh);
[Vx,Vy,~,~,~,~]=holomorphic_cont_tri(Tinv(:,equationnumber),holo_arc,Ppm,Qqm,Lf,Uf,Start,I,degree,bus_n);
% [Vx,Vy,~,~,~,~]=holomorphic_cont_tri_mex(Tinv(:,equationnumber),holo_arc,Ppm,Qqm,Lf,Uf,Start,I,degree,bus_n);

[num,den]=Pade_Apprxmt(Vx,Vy,degree,bus_n);
% [num,den]=Pade_Apprxmt_mex(Vx,Vy,degree,bus_n);

[~,~,xholo,dba]=update_V_Pade(num,den,degree,Tinv(:,equationnumber),holo_arc,arc_ratio,bus_n,I);
% [~,~,xholo,dba]=update_V_Pade_mex(num,den,degree,Tinv(:,equationnumber),holo_arc,arc_ratio,bus_n,I);

holo_err=quad_form_vals(xholo,Ma,2*bus_n-1)-bat-dba;

if aal_ind==0
   % Pade pole estimate by CONSENSUS across a fixed set of buses. The true singularity
   % sits at ~the same magnitude on every bus; spurious (Froissart) poles are bus-specific.
   % Pick the smallest pole magnitude agreed on by >= half the sampled buses; if none
   % reaches consensus, fall back to the global smallest (original behaviour).
   sampled_buses=param.pole_buses(param.pole_buses>=1 & param.pole_buses<=bus_n);
   nb=numel(sampled_buses);
   pmag=zeros(0,1); pbus=zeros(0,1);
   for ii=1:nb
       r=roots(flip(den(sampled_buses(ii),:)));
       v=isfinite(r) & sign(real(r))==sgn & abs(imag(r))<param.im_tol;
       rv=abs(r(v));
       pmag=[pmag; rv]; pbus=[pbus; ii*ones(numel(rv),1)]; %#ok<AGROW>
   end
   if ~isempty(pmag)
       [sm,ord]=sort(pmag); sb=pbus(ord); N=numel(sm);
       if param.pole_consensus
           minsupp=max(2,round(nb/2)); sel=NaN; i=1;
           while i<=N
               j=i;
               while j<N && (sm(j+1)-sm(i))<=param.pole_cl_tol*sm(i); j=j+1; end
               if numel(unique(sb(i:j)))>=minsupp; sel=median(sm(i:j)); break; end
               i=j+1;
           end
           if isnan(sel); sel=sm(1); end
       else
           sel=sm(1);
       end
       aal_limit=aal+sgn*sel*holo_arc;
       aal_ind=1;
   end
end

max_mismatch=max(abs(holo_err));
kdst1=0;
while max_mismatch<param.mismatch_error && sgn*(aal+arc_ratio*holo_arc)<=sgn*aal_limit
    kdst1=kdst1+1;
    arc_ratio=arc_ratio*param.holo_arc_grow;
    [~,~,xholo,dba]=update_V_Pade(num,den,degree,Tinv(:,equationnumber),holo_arc,arc_ratio,bus_n,I);
%     [~,~,xholo,dba]=update_V_Pade_mex(num,den,degree,Tinv(:,equationnumber),holo_arc,arc_ratio,bus_n,I);
    holo_err=quad_form_vals(xholo,Ma,2*bus_n-1)-bat-dba;
    max_mismatch=max(abs(holo_err));
    if Display==1
        fprintf('Increase arc length! Step %d, max_mismatch = %g\n',kdst1,max_mismatch);
    end
    if abs(arc_ratio*holo_arc)>param.holo_unbounded
        fprintf('Unbounded trace! Skip. \n');
        break;
%     elseif abs(arc_ratio*holo_arc)>=6.*10^(-3)
%         fprintf('Reach maximum arc length! \n');
%         arc_ratio=sgn*6*10^(-3)/holo_arc;
%         break;
    end
end

kdst2=0;
while max_mismatch>=param.mismatch_error || sgn*(aal+arc_ratio*holo_arc)>=sgn*aal_limit
    kdst2=kdst2+1;
    arc_ratio=arc_ratio*param.holo_arc_shrink;
    [~,~,xholo,dba]=update_V_Pade(num,den,degree,Tinv(:,equationnumber),holo_arc,arc_ratio,bus_n,I);
%     [~,~,xholo,dba]=update_V_Pade_mex(num,den,degree,Tinv(:,equationnumber),holo_arc,arc_ratio,bus_n,I);
    holo_err=quad_form_vals(xholo,Ma,2*bus_n-1)-bat-dba;
    max_mismatch=max(abs(holo_err));
    if Display==1
        fprintf('Decrease arc length! Step %d, max_mismatch = %g\n',kdst2,max_mismatch);
    end
    if abs(arc_ratio*holo_arc)<param.aal_min || holo_arc/holo_arcp<=param.holo_switch_ratio
        flagh=1;
        if Display==1
            fprintf('Too small holomorphic convergence range! Switch to numerical continuation. \n');
        end
        break;
    end
end
holo_arc1=abs(arc_ratio*holo_arc);
aal1=aal+arc_ratio*holo_arc;

end
