function [Z,Y,Fail,sgn,cmplt]=branch_trace_hybrid_4_no_mex(x,x0,M0,MS,km0,bbt,bb,...
    svec,param,sgn0,aal,Dspl,Fail,Y_start,purpose,Tinv,escape,solutionnumber,equationnumber)

% clear F J Y Z

sarc = zeros(1,param.smax+3);  % arclength.
Fail(solutionnumber,equationnumber)=-1;
cross_singular=0; % index for across singularity
% countdown=param.countdown;
cmplt=0;
% param.interval=5;
% param.slope_max=2*10^2;

Z = []; 
z = zeros(size(x));
Nm = size(x,1);
Nx = size(x,1);
Y=zeros(Nx+1,param.smax+3);
fid1=0;
fid2=0;
    
if purpose==1 % cold start
    alpha = param.alpha_init; % initial alpha
    Y(:,1) = [x; 0];  % initial state;
    xm1 = x; 
    am1 = alpha;
%% first step
    for icount = 1:param.first_step_imax
        b=2*M0*x;
        d0=x'*b/2;
        % D variant: Jacobian built directly via MakeJacobianD below
        J=MakeJacobianD(x,MS,Nm);
        F = J*x/2+d0*km0-bbt - svec*alpha;
    
        J1=[J param.fac*km0;b'/param.fac^2 -1/param.fac];
        F1=[F;0];
        dxx=-J1\F1;
        dx=dxx(1:Nx);
        x = x+dx;
        if Dspl==1
            fprintf('Step 1, |dx| = %g, |F| = %g\n',norm(dx),norm(F));
        end
        if (norm(F)<param.ef)&&(norm(dx,inf)<param.ex)
            break;
        end
    end

    if icount>=param.first_step_imax
        alpha = -alpha/param.alpha_retry_div; % initial alpha
        x=Y(1:end-1,1);  % initial state;
        xm1 = x; 
        am1 = alpha;
        % first step second trial
        for icount2 = 1:param.first_step_imax
            b=2*M0*x;
            d0=x'*b/2;
            % D variant: Jacobian built directly via MakeJacobianD below
            J=MakeJacobianD(x,MS,Nm);
            F = J*x/2+d0*km0-bbt - svec*alpha;
    
            J1=[J param.fac*km0;b'/param.fac^2 -1/param.fac];
            F1=[F;0];
            dxx=-J1\F1;
            dx=dxx(1:Nx);
            x = x+dx;
            if Dspl==1
                fprintf('Step 1, |dx| = %g, |F| = %g\n',norm(dx),norm(F));
            end
            if (norm(F)<param.ef)&&(norm(dx,inf)<param.ex)
                break;
            end
        end
        
        if icount2==param.first_step_imax
            fid1=1;
            fprintf('\nFail to obtain the first step! Stop!\n')
        end
    end
%% second step
    if fid1==0 % if the first step is solved, try the second step
        ds = [(x-xm1); alpha];  % initial step
        Y(:,2) = [x; alpha];
        snew = sarc(1) + norm(ds);
        sarc(2) = snew;

% set  next step size
        astep = param.sstep/norm(ds);
        astep = min(astep,param.maxphase2);  % set a maximum stepsize
        astep = max(astep,param.minimumstep/param.second_step_minstep_div);
% astep = max(astep,param.minimumstep);

        xm1 = x; 
        am1 = alpha;
% predictor
        x = x + ds(1:Nm)*astep; 
        alpha = alpha + ds(end)*astep;
% corrector
        for icount = 1:param.first_step_imax
            b=2*M0*x;
            d0=x'*b/2;
            % D variant: Jacobian built directly via MakeJacobianD below
            J=MakeJacobianD(x,MS,Nm);
            F = J(:,1:Nm)*x/2+d0*km0-bbt - svec*alpha;
            J(:,Nm+1) = -svec;
    
            ds = [(x-xm1); (alpha-am1)];   
            F(Nm+1,1) = (ds'*ds) - param.sstep^2;
            J(Nm+1,:) = [2*(x-xm1)' 2*(alpha-am1)];   
    
            J1=sparse([J [param.fac*km0;0];[b'/param.fac^2 0] -1/param.fac]);
            F1=[F;0];
            dy=-J1\F1;
            dx=dy(1:Nx);
            x = x+dx; 
            alpha = alpha + dy(Nx+1);
            if Dspl==1
                fprintf('Step 2, |dx| = %g, |F| = %g, Stepsize = %g\n',norm(dx,inf),norm(F),astep);
            end
            if (norm(F)<param.ef)&&(norm(dx,inf)<param.ex)
                break;
            end
        end

        if icount>=param.first_step_imax
            fid2=1;
            fprintf('\nFail to obtain the second step! Stop!\n')
        end
% icount
        ds = [(x-xm1); (alpha-am1)];  % actual step
        dsm1 = ds;
        Y(:,3) = [x; alpha];
        snew = sarc(2) + norm(ds);
        sarc(3) = snew;

        actstep=3;
        scount = 1;
% astep=astep+astepp*0.03;
    end
elseif purpose==0
    actstep=3;
    scount = 1;
    Y(:,1:3)=Y_start;
    x=Y(1:end-1,3);
    sarc(2)=norm(Y(:,2)-Y(:,1))+sarc(1);
    sarc(3)=norm(Y(:,3)-Y(:,2))+sarc(2);
    alpha=Y(end,3);
    astep = norm(Y(:,3)-Y(:,2))/param.sstep;
    astep = min(astep,param.maxphase1);  % maximum step
    if isfield(param,'init_astep_cap'); astep=min(astep,param.init_astep_cap); end % fold-proximity cap on hand-back initial step
    astep = max(astep,param.minimumstep); % minstep
end

%% continuation
% Adaptive per-step cap (param.max_adapt): adaptive maximumstep = param.maximumstep * ad_ratio,
% ad_ratio = slope_prev/slope (secant-slope change) capped at 1, recomputed FRESH each step on
% the fixed base (NOT compounded) so the ceiling snaps back to full when the slope steadies.
% slope = 2-point secant |Dstate|/Dalpha. Shrinks the step approaching a fold, restores it after.
slope_prev = [];
if (purpose==1 && fid1+fid2==0) || purpose==0

for outersteps = 1:param.smax
    
%     if cross_singular==1
%         countdown=countdown-1;
%     end
%     
%     if countdown<=0
%         break;
%     end
    
    xm1 = x; 
    am1 = alpha;
    sm0 = sarc(actstep);    ym0 = Y(:,actstep);
    sm1 = sarc(actstep-1);  ym1 = Y(:,actstep-1);
    sm2 = sarc(actstep-2);  ym2 = Y(:,actstep-2);
    dsm1=sm1-sm0;
    dsm2=sm2-sm0;
    sp1 = sm0 + astep*param.sstep;

    spars1=ym0;
    spars2 = (dsm1^2*(ym2-ym0)-dsm2^2*(ym1-ym0))/(dsm1*dsm2*(dsm1-dsm2));
    spars3 = (dsm2*(ym1-ym0)-dsm1*(ym2-ym0))/(dsm1*dsm2*(dsm1-dsm2));
    x = spars1(1:Nx) + spars2(1:Nx)*(sp1-sm0) + spars3(1:Nx)*(sp1-sm0)^2;
    alpha = spars1(1+Nx) + spars2(1+Nx)*(sp1-sm0) + spars3(1+Nx)*(sp1-sm0)^2;
    actstep=actstep+1;  
%}   
    
    % corrector
    for icount = 1:param.imax
        b=2*M0*x;
        d0=x'*b/2;
        % D variant: Jacobian built directly via MakeJacobianD below
        J=MakeJacobianD(x,MS,Nm);
        F = J(:,1:Nm)*x/2+d0*km0-bbt - svec*alpha;
        J(:,Nm+1) = -svec;

        ds = [(x-xm1); (alpha-am1)];
        F(Nm+1,1) = (ds'*ds) - (astep*param.sstep)^2;
        J(Nm+1,:) = 2*ds';

        J1=[J [param.fac*km0;0];[b'/param.fac^2 0] -1/param.fac];
        F1=[F;0];
        if param.Dind == 0 % phase I
           dy=-J1\F1;
           dx=dy(1:Nx);
           x = x+ dx;
           alpha = alpha + dy(Nx+1);
           if (norm(F)<param.ef)&&(norm(dx,inf)<param.ex)
              break;
           end
        else
            if param.Dind ~= 0 % phase II
               param.maximumstep = param.maxphase2;
               PhaseStep=PhaseStep-1;
               if PhaseStep == 0 % return to phase I
                  param.Dind=0; 
                  param.maximumstep = param.maxphase1;
%                   fprintf('Return to Phase I.\n')
%                   fprintf('=================================================================\n\n')
               end
            end
            J1=full(J1);
            dy=-pinv(J1)*F1;
           dx=dy(1:Nx);
           x = x+ dx;
           alpha = alpha + dy(Nx+1);
           if (norm(F)<param.phase2_ef_mult*param.ef)&&(norm(dx,inf)<param.phase2_ex_mult*param.ex)
              break;
           end
        end
    end  
    if Dspl==1
        rmd=rem(outersteps,param.dev);% interval display
        if rmd==1
            fprintf('Trace # (%d,%d), Skip_ind=%d, a=%g, Step # %d, Stepsize: %g, # Iter: %d, # Solu: %d\n',solutionnumber,equationnumber,param.skipindex,alpha,actstep,astep,icount,scount);
        end
    end
    %}
     if icount > param.stall_iter_mult*param.targetsteps
        param.skipindex=param.skipindex+1; 
        if param.skipindex > param.maxskip % reach the maximum number of steps
            Fail(solutionnumber,equationnumber)=1;
%             fprintf('===================================================\n');
%             fprintf('Too many numerical stalls! Skip! \n');
%             fprintf('===================================================\n'); 
            break;
        end
        if astep<=param.phase2_trigger_mult*param.minimumstep
            param.Dind=1;
            PhaseStep=param.PS;
%             fprintf('\n=================================================================\n')
%             fprintf('Numerically unstable! Back up and execute Phase II for %g steps.\n',param.PS)
        else
%             fprintf('\n=================================================================\n')
%             fprintf('Numerically unstable! Back up.\n')
        end
        % back up several steps, decrease step
        ds = [(x-xm1); (alpha-am1)];  % actual step
        Y(:,actstep)=[x; alpha];
        snew = sarc(actstep-1) + norm(ds);
        sarc(actstep) = snew;
        dsm1 = ds;
         astep = astep*(1 + (param.targetsteps-icount)/param.imax);
         astep = min(astep,param.maximumstep);  % maximum step
         astep = max(astep,param.minimumstep); % minstep
         
        if actstep==4
            if icount<param.backup_iter_frac*param.imax
                actstep=actstep-1;
                sarc(actstep+1:actstep+1) = 0;
                Y(:,actstep+1:actstep+1)=zeros(Nx+1,1);
                x = Y(1:end-1,actstep); 
                alpha = Y(end,actstep);
                astep = astep/param.backup_div_small;
            else
                actstep=actstep-1;
                sarc(actstep+1:actstep+1) = 0;
                Y(:,actstep+1:actstep+1)=zeros(Nx+1,1);
                x = Y(1:end-1,actstep); 
                alpha = Y(end,actstep);
                astep = astep/param.backup_div_large;
            end          
            astep = max(astep,param.minimumstep);
        elseif actstep>4
            if icount<param.backup_iter_frac*param.imax
                actstep=actstep-1;
                sarc(actstep+1:actstep+1) = 0;
                Y(:,actstep+1:actstep+1)=zeros(Nx+1,1);
                x = Y(1:end-1,actstep); 
                alpha = Y(end,actstep);
                astep = astep/param.backup_div_small;
            else
                actstep=actstep-2;
                sarc(actstep+1:actstep+2) = 0;
                Y(:,actstep+1:actstep+2)=zeros(Nx+1,2);
                x = Y(1:end-1,actstep); 
                alpha = Y(end,actstep);
                astep = astep/param.backup_div_large;
            end          
            astep = max(astep,param.minimumstep);
        end       
    else
        ds = [(x-xm1); (alpha-am1)];  % actual step
        Y(:,actstep)=[x; alpha];
        snew = sarc(actstep-1) + norm(ds);
        sarc(actstep) = snew;
        dsm1 = ds;
         astep = astep*(1 + (param.targetsteps-icount)/param.step_adapt_denom);
         astep = min(astep,param.maximumstep);  % maximum step
         if isfield(param,'max_adapt') && param.max_adapt   % adaptive cap = base * ad_ratio(secant-slope change)
             dstate=norm(x-xm1); dalpha=abs(alpha-am1);
             if dalpha>0
                 slope=dstate/dalpha;
                 if ~isempty(slope_prev) && slope>0
                     ad_ratio = min(slope_prev/slope, 1);
                     astep = min(astep, max(param.maximumstep*ad_ratio, param.minimumstep));
                 end
                 slope_prev=slope;
             end
         end
         astep = max(astep,param.minimumstep); % minstep
        if (aal+sgn0*alpha)*(aal+sgn0*am1)<=0  % near solution.  find it.
%             fprintf('========================================\n');
%             fprintf('Near solution. Find it.\n');
            z = abs(alpha)/(abs(alpha)+abs(am1))*xm1 + abs(am1)/(abs(alpha)+abs(am1))*x;
            scount = scount+1;
            for zcount = 1:param.imax
                bz=2*M0*z;
                d0=z'*bz/2;
                
                Jz=MakeJacobianD(z,MS,Nm);
                Fz = Jz*z/2+d0*km0-bb;
                
                Jz1=full([Jz param.fac*km0;bz'/param.fac^2 -1/param.fac]);
                Fz1=[Fz;0];
                dzz = -pinv(Jz1)*Fz1;
                dz=dzz(1:Nx);
                z = z+dz;
                if Dspl==1
                    fprintf('|dz| = %g, |Fz| = %g\n',norm(dz),norm(Fz));
                end
                if (norm(Fz)<param.ef) && (norm(dz,inf)<param.ex)
                    Z = [Z z];
%                     fprintf('P-C finds a solution!\n');
%                     fprintf('========================================\n');
                    break;
                end
                
                if zcount==param.imax
%                     fprintf('Fail to resolve the solution! Back up!\n');
%                     fprintf('========================================\n');
                    actstep=actstep-1;
                    sarc(actstep+1:actstep+1) = 0;
                    Y(:,actstep+1:actstep+1)=zeros(Nx+1,1);
                    x = Y(1:end-1,actstep); 
                    alpha = Y(end,actstep);
                    astep = astep/param.zcount_fail_div;
                end
            end
%             Z = [Z z];
            
            if scount>param.maxsolu
                Fail(solutionnumber,equationnumber)=2; % number of solutions over maximum number
%                 fprintf('===================================================\n');
%                 fprintf('Too many solutions at this loop! Skip and continue!\n');
%                 fprintf('===================================================\n'); 
                break;
            elseif (norm(z-x0,inf)<param.ex)
%                 fprintf('========================================================\n');
%                 fprintf('Trace completed!\n');
%                 fprintf('========================================================\n');
                Fail(solutionnumber,equationnumber)=0; % succeed!
                cmplt=1;
                break;
            end            
        end       
     end
     
     if actstep>=1+2*param.interval && (Y(end,actstep)-Y(end,actstep-param.interval))*(Y(end,actstep-param.interval)-Y(end,actstep-2*param.interval))<0
         cross_singular=1;
         paal=Y(end,actstep-param.interval);
%          sgn=sign(Y(end,actstep)-Y(end,actstep-10));       
     end
     
     if cross_singular==1
         secant_slope=norm((Y(1:end-1,actstep)-Y(1:end-1,actstep-1))/(Y(end,actstep)-Y(end,actstep-1)),inf);
         if secant_slope<=param.slope_max && abs(Y(end,actstep)-paal)>=param.handback*param.aal_min % turning-point hand-back hysteresis (param.handback): clears the near-singular zone before returning to holomorphic, so the curve doesn't get misled onto a wrong/non-closing branch
             break;
         end
     end
     
     if escape==1 % escape from zero constant point
         dba=Tinv*(bbt+svec*alpha);
         if min(abs(dba))>=param.escape_tol % 1.5*10^(-4)
             secant_slope=norm((Y(1:end-1,actstep)-Y(1:end-1,actstep-1))/(Y(end,actstep)-Y(end,actstep-1)),inf);
             if secant_slope<=param.slope_max
                break;
             end
         end
     end
end

Y = Y(:,1:actstep); % actual non-zero values
Y(end,:)=Y(end,:)*sgn0;
sgn=sign(Y(end,actstep)-Y(end,actstep-1));

if outersteps == param.smax % reach the maximum number of steps
   Fail(solutionnumber,equationnumber)=3;
%    fprintf('===================================================\n');
%    fprintf('Reach iteration max number! \n');
%    fprintf('===================================================\n'); 
end
%
[~,Kh]=size(Z);
rpt=0;

Dt=[];
for i=1:Kh-1
     for j=i+1:Kh
         if norm(Z(:,j))~=0 && ( norm(Z(:,j)-Z(:,i),inf)<param.ex || norm(Z(:,j)+Z(:,i),inf)<param.ex )
             Z(:,j)=0;
             rpt=rpt+1;
             Dt(1,rpt)=j;
         end
     end
end

Z(:,Dt)=[];
%}
end