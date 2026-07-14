%% trace loops
% clear F J Y Z

% global fac bus_n numofcons

%% general settings
param.sstep = 1e-5;
param.targetsteps = 3;                  % target number of iterations.
param.imax = 40;                        % maximum number of Newton iterations
param.smax = 20000;                     % maximum number of steps
param.maximumstep = 8000;               % maximum step length
param.maxphase1 = param.maximumstep/2;  % maximum step length for phase I
param.maxphase2 = param.maxphase1/10;   % maximum step length for phase II
param.minimumstep = 1e-2;               % minimum step length
param.dev=20;                            % denominator of interval printing
param.skipindex=0;                      % skip index
param.ef=1e-8;                          % tolerence for function value
param.ex=1e-7;                          % tolerence for dx value
param.Dind=0;                           % trigger for changing phases
param.PS=300;                           % phase II total steps
param.condn=1e-8;                       % condition number
param.maxsolu=100;                      % maximum number of solutions for each path
param.maxskip=50;                       % maximum number of numerical stalling
% param.countdown=10;                     % after crossing singularity how far it will go with numerical continuation
param.fac=fac;
warning off MATLAB:nearlySingularMatrix;

% Validate and initialize case-specific sparse operators once per trace.
MakeJacobianD([],MS,numofcons,'initialize');
quad_form_vals([],Ma,numofcons,'initialize');

% Make the random two-bus Pade sampling reproducible and independent of
% parallel worker scheduling or solution discovery order.
seed_weights=(1:numel(x))';
solution_signature=sum(abs(round(x*1e8)).*seed_weights);
trace_seed=mod(solution_signature+104729*bus_n+1009*equationnumber,2^32-1);
rng(double(trace_seed),'twister');

%% special settings
param.mismatch_error=0.1*10^-3; % deafault 0.1*10^-3
param.holonum=50;               % default 50
param.holo_arc_max=0.01;        % default 0.02
param.outnum=500;               % default 500
param.slope_max=4*10^5;         % shared default (4e5)
param.interval=1;               % default 1
param.aal_min=0.1*10^(-5);      % default 0.5*10^(-5)
param.im_tol=1*10^(-3);         % default 1*10^(-3)
param.handback = 0.1; % It affects the numerical stability in different systems and curves. For case3 and case118, 0.5 is a better choice. For case39 and case30, 0.1 is a good choice.

%% consensus pole detection (see holomorphic_hybrid_5_no_mex.m)
% NOTE: the per-curve RNG is already seeded above (trace_seed / rng), so the 'random'
% bus set is REPRODUCIBLE per (solution,equation) and independent of worker scheduling.
global POLEMODE %#ok<GVMIS>  % experiment override: 'fixed' | 'random' (default 'fixed')
if ~isempty(POLEMODE); param.pole_mode=POLEMODE; else; param.pole_mode='fixed'; end
param.pole_fixed_buses = 1:2;   % option 1: predefined bus set (auto-filtered to <=bus_n)
param.pole_nbus        = 5;     % option 2: size of the random bus set
param.pole_consensus   = true;  % pick smallest pole agreed on by >=half the buses (else global-min)
param.pole_cl_tol      = 0.10;  % relative magnitude tolerance for clustering "the same pole"
% resolve the bus set ONCE PER CURVE (held fixed for the whole trace)
if strcmp(param.pole_mode,'random')
    param.pole_buses = sort(randperm(bus_n, min(param.pole_nbus, bus_n)));
else
    param.pole_buses = param.pole_fixed_buses(param.pole_fixed_buses<=bus_n & param.pole_fixed_buses>=1);
    if isempty(param.pole_buses); param.pole_buses = 1:min(param.pole_nbus,bus_n); end
end

%% display settings
holo_display=0;
reso_display=0;
cont_display=0;

%% initialization
%
param.nsolu=solutionnumber;
param.neqt=equationnumber;
%}
%{
param.nsolu=13;
param.neqt=2;
svec=zeros(numofcons,1); 
svec(param.neqt) = 1;
Tholo=0;
Tresol=0;
Tcont=0;
Tdata=0;
NYh=0;
Fail=[];
x=Zsave(:,param.nsolu);
%} 
purpose=0;                              % 0: warm start, 1: cold start
nump=zeros(bus_n,degree.numerator);
denp=zeros(bus_n,degree.denominator);
arc_ratiop=0;
Y=[x;0];
Z=x;
sgn=1;
bat=ba;
cmplt=0;
cmplt1=0;
escape=0;                               % escape from 0 constant

%% if need to escape from 0 constant
if min(abs(ba))<10^(-6)
    purpose=1;
    escape=1;
    xstart=x;
    aalstart=0;
    bbt=bb;
    Y_start=zeros(numofcons+1,3);
    
    [Z1,Y1,Fail,sgn,cmplt1]=...
        branch_trace_hybrid_4_no_mex(xstart,x,M0,MS,km0,bbt,bb,sgn*svec,...
        param,sgn,aalstart,cont_display,Fail,Y_start,purpose,Tinv,escape,solutionnumber,equationnumber); 
    Z=[Z Z1];
    Y = [Y Y1];
end

%% main body
for outstep=1:param.outnum
    if cmplt1==1                        % check if the trace has been completed
       break; 
    end
    
    if outstep==1 && escape==0          % initialize the first step without 0 constant
        holo_arc=zeros(param.holonum,1);
        holo_arc(1)=param.holo_arc_max/2;
        aal=zeros(param.holonum,1);
        xholo=zeros(numofcons,param.holonum);
        xholo(:,1)=x;
    elseif outstep==1 && escape==1      % initialize the first step with 0 constant
        holo_arc=zeros(param.holonum,1);
        holo_arc(1)=param.holo_arc_max/2;
        aal=zeros(param.holonum,1);
        aal(1)=Y(end,end);
        xholo=zeros(numofcons,param.holonum);
        xholo(:,1)=Y(1:end-1,end);
        escape=0;
    else                                % initialize the next step from the last step
        holo_arc=zeros(param.holonum,1);
        holo_arc(1)=hal;
        aal=zeros(param.holonum,1);
        aal(1)=Y(end,end);
        xholo=zeros(numofcons,param.holonum);
        xholo(:,1)=Y(1:end-1,end);
    end   
    Y_start=zeros(numofcons+1,3);    
    aal_limit=100*sgn;                  % iniliaze the nearest pole distance
    aal_ind=0;                          % iniliaze the index which indicates when to stop computing poles
    
    for holostep=1:param.holonum        % start holomorphic step
        if holo_display==1              % display holomorphic steps
            fprintf('Outer step: %d, Holo step: %d\n',outstep,holostep);
        end
        txy=xholo(bus_n+1:end,holostep);
        txy1=zeros(bus_n,1);
        txy1(1:I.slack-1)=txy(1:I.slack-1);
        txy1(I.slack+1:end)=txy(I.slack:end);
        solu0=[xholo(1:bus_n,holostep);txy1];
        dba=Tinv(:,param.neqt)*aal(holostep);
        bat=ba+dba;
    %% find holomorphic prediction
        if holostep==1 && outstep==1    % initialize the first holomorphic step at the first outer step
            ha=holo_arc(holostep);
            hap=param.aal_min;
        elseif holostep==1 && outstep>1 % initialize the first holomorphic step beyond the first outer step
            ha=50*param.aal_min;
            hap=param.aal_min;
        elseif holostep==2              % initialize the second holomorphic step
            ha=holo_arc(holostep-1);
            hap=param.aal_min;
        else                            % initialize the rest holomorphic step
            if holo_arc(holostep-1)>1.5*holo_arc(holostep-2)
                ha = holo_arc(holostep-1);
            else
                ha=holo_arc(holostep-1)/2;
            end
            hap=holo_arc(holostep-1);
        end
        
        [xnew,aal(holostep+1),aal_limit,aal_ind,holo_arc(holostep),arc_ratio,num,den,kdst,flagh]=...
            holomorphic_hybrid_5_no_mex(Mp,Mq,solu0,Ybus,I,Ma,bat,Tinv,ha,hap,degree,...
            param,aal(holostep),aal_limit,aal_ind,sgn,bus_n,equationnumber,holo_display);
        
        %% holomorphic progress small, switch to numerical continuation
        if flagh==1                     % if holomorphic progress is too small
            if holostep>1               % if the holomorphic step is beyond the first step, initialize the predictor-corrector algorithm by a warm starter
                purpose=0;
                xholo(:,holostep+1:end)=[];
                holo_arc(holostep:end)=[];
                aal(holostep+1:end)=[];
                if abs(aal_limit)<100
                    aal_interval=abs(aal_limit-aal(holostep));
                    krt=(aal_interval/5)/holo_arc(holostep-1);
                    if krt>0.45
                        krt=0.45;
                    end
                else
                    krt=0.05;
                end
                bbt=bb;
                bbt(param.neqt)=bbt(param.neqt)+aal(holostep-1)+sgn*(1-2*krt)*holo_arc(holostep-1);
                [~,~,Y_start(1:end-1,1),~]=update_V_Pade(nump,denp,degree,Tinv(:,param.neqt),holo_arc(holostep-1),(1-2*krt)*arc_ratiop,bus_n,I);
                [~,~,Y_start(1:end-1,2),~]=update_V_Pade(nump,denp,degree,Tinv(:,param.neqt),holo_arc(holostep-1),(1-krt)*arc_ratiop,bus_n,I);
%                 [~,~,Y_start(1:end-1,1),~]=update_V_Pade_mex(nump,denp,degree,Tinv(:,param.neqt),holo_arc(holostep-1),(1-2*krt)*arc_ratiop,bus_n,I);
%                 [~,~,Y_start(1:end-1,2),~]=update_V_Pade_mex(nump,denp,degree,Tinv(:,param.neqt),holo_arc(holostep-1),(1-krt)*arc_ratiop,bus_n,I);
                Y_start(1:end-1,3)=xholo(:,holostep);
                Y_start(end,1)=0;
                Y_start(end,2)=krt*holo_arc(holostep-1);
                Y_start(end,3)=2*krt*holo_arc(holostep-1);                          
                xstart=xholo(:,end);
                aalstart=aal(end-1)+sgn*(1-2*krt)*holo_arc(end);
            else                        % if the holomorphic step is the first step, initialize the predictor-corrector algorithm by a cold starter
                purpose=1;
                xholo(:,holostep+1:end)=[];
                holo_arc(holostep+1:end)=[];
                holo_arc(holostep)=3*param.aal_min;
                aal(holostep+1:end)=[];               
                xstart=xholo;
                aalstart=aal;
                bbt=bb;
                bbt(param.neqt)=bbt(param.neqt)+aal;
            end
            break;
        else                            % if holomorphic progress is good, continue it
            nump=num;                   % record step Pade numerator
            denp=den;                   % record step Pade denominator
            arc_ratiop=arc_ratio;       % record step arc_ratio
        end
        %% correct holomorphic prediction
        bbt=bb;
        bbt(param.neqt)=bbt(param.neqt)+aal(holostep+1);
        
        [xnew,flagr] = Resolve_no_mex(MS,bbt,M0,km0,xnew,numofcons,reso_display);
        
        if flagr==1                     % success of correcting
            xholo(:,holostep+1)=xnew;   % recoganize the new step as valid
            numpp=num;                  % record step Pade numerator
            denpp=den;                  % record step Pade denominator
            arc_ratiopp=arc_ratio;      % record step arc_ratio
        end
        
        %% jump out criteria
        if flagr==-1                    % correction failure, delete current point and switch to numerical continuation
            if holostep<2               % if holomorphic step is the first step, initialize predictor-corrector algorithm from a cold starter
                purpose=1;
                xholo(:,holostep+1:end)=[];
                holo_arc(holostep:end)=[];
                aal(holostep+1:end)=[];
                bbt=bb;
                bbt(param.neqt)=bbt(param.neqt)+aal(holostep);                
                xstart=xholo(:,end);
                aalstart=aal(end);
            else                        % if holomorphic step is not the first step, initialize predictor-corrector algorithm from a warm starter
                purpose=0;
                xholo(:,holostep+1:end)=[];
                holo_arc(holostep:end)=[];
                aal(holostep+1:end)=[];
                if abs(aal_limit)<100
                    aal_interval=abs(aal_limit-aal(holostep));
                    krt=(aal_interval/5)/holo_arc(holostep-1);
                    if krt>0.45
                        krt=0.45;
                    end
                else
                    krt=0.03;
                end
                bbt=bb;
                bbt(param.neqt)=bbt(param.neqt)+aal(holostep-1)+sgn*(1-2*krt)*holo_arc(holostep-1);
                [~,~,Y_start(1:end-1,1),~]=update_V_Pade(numpp,denpp,degree,Tinv(:,param.neqt),holo_arc(holostep-1),(1-2*krt)*arc_ratiopp,bus_n,I);
                [~,~,Y_start(1:end-1,2),~]=update_V_Pade(numpp,denpp,degree,Tinv(:,param.neqt),holo_arc(holostep-1),(1-krt)*arc_ratiopp,bus_n,I);
%                 [~,~,Y_start(1:end-1,1),~]=update_V_Pade_mex(numpp,denpp,degree,Tinv(:,param.neqt),holo_arc(holostep-1),(1-2*krt)*arc_ratiopp,bus_n,I);
%                 [~,~,Y_start(1:end-1,2),~]=update_V_Pade_mex(numpp,denpp,degree,Tinv(:,param.neqt),holo_arc(holostep-1),(1-krt)*arc_ratiopp,bus_n,I);
                Y_start(1:end-1,3)=xholo(:,holostep);
                Y_start(end,1)=0;
                Y_start(end,2)=krt*holo_arc(holostep-1);
                Y_start(end,3)=2*krt*holo_arc(holostep-1);                          
                xstart=xholo(:,end);
                aalstart=aal(end-1)+sgn*(1-2*krt)*holo_arc(end);
            end
            %}
            break;
        end
        
        %% near solution, find it
        if aal(holostep+1)*aal(holostep)<0  % identify real solution
%             fprintf('========================================\n');
%             fprintf('Cross solution. Find it.\n');
            arc_ratio1=arc_ratiop*abs(aal(holostep)/(aal(holostep)-aal(holostep+1)));
            
            [~,~,z,~]=update_V_Pade(nump,denp,degree,Tinv(:,param.neqt),holo_arc(holostep),arc_ratio1,bus_n,I);
%             [~,~,z,~]=update_V_Pade_mex(nump,denp,degree,Tinv(:,param.neqt),holo_arc(holostep),arc_ratio1,bus_n,I);

            [znew,flagr2] = Resolve_no_mex(MS,bb,M0,km0,z,numofcons,reso_display);
            
            if norm(znew-x,inf)<10^(-7)     % check completeness of the trace
                cmplt=1;
                break;
            end
            
            if flagr2==1                    % if resolving solution succeeds, record new solution
                Z=[Z znew];
%                 fprintf('HE finds a solution!\n');
%                 fprintf('========================================\n');
            else                            % if resolving solution fails, indicating that last holomorphic prediction is wrong. Back up and run predictor-corrector algorithm from a cold starter
                purpose=1;
                xholo(:,holostep+1:end)=[];
                holo_arc(holostep:end)=[];
                aal(holostep+1:end)=[];
                bbt=bb;
                bbt(param.neqt)=bbt(param.neqt)+aal(holostep);                
                xstart=xholo(:,end);
                aalstart=aal(end);
                %{
                if holostep<2
                    purpose=1;
                    xholo(:,holostep+1:end)=[];
                    holo_arc(holostep:end)=[];
                    aal(holostep+1:end)=[];
                    bbt=bb;
                    bbt(param.neqt)=bbt(param.neqt)+aal(holostep);                
                    xstart=xholo(:,end);
                    aalstart=aal(end);
                else
                    purpose=0;
                    xholo(:,holostep+1:end)=[];
                    holo_arc(holostep:end)=[];
                    aal(holostep+1:end)=[];
%                     if abs(aal_limit)<100
%                         aal_interval=abs(aal_limit-aal(holostep));
%                         krt=(aal_interval/5)/holo_arc(holostep-1);
%                         if krt>0.45
%                             krt=0.45;
%                         end
%                     else
                        krt=0.01;
%                     end
                    bbt=bb;
                    bbt(param.neqt)=bbt(param.neqt)+aal(holostep-1)+sgn*(1-2*krt)*holo_arc(holostep-1);
%                 [~,~,Y_start(1:end-1,1),~]=update_V_Pade(numpp,denpp,degree,Tinv(:,param.neqt),holo_arc(holostep-1),(1-2*krt)*arc_ratiopp,bus_n,I);
%                 [~,~,Y_start(1:end-1,2),~]=update_V_Pade(numpp,denpp,degree,Tinv(:,param.neqt),holo_arc(holostep-1),(1-krt)*arc_ratiopp,bus_n,I);
                    [~,~,Y_start(1:end-1,1),~]=update_V_Pade(numpp,denpp,degree,Tinv(:,param.neqt),holo_arc(holostep-1),(1-2*krt)*arc_ratiopp,bus_n,I);
                    [~,~,Y_start(1:end-1,2),~]=update_V_Pade(numpp,denpp,degree,Tinv(:,param.neqt),holo_arc(holostep-1),(1-krt)*arc_ratiopp,bus_n,I);
                    Y_start(1:end-1,3)=xholo(:,holostep);
                    Y_start(end,1)=0;
                    Y_start(end,2)=krt*holo_arc(holostep-1);
                    Y_start(end,3)=2*krt*holo_arc(holostep-1);                          
                    xstart=xholo(:,end);
                    aalstart=aal(end-1)+sgn*(1-2*krt)*holo_arc(end);
                end
                %}
                break;
            end
        end
    end
    
    if isempty(holo_arc)        % check if there is a successful holomorphic step, if not
        hal=param.holo_arc_max/20;
    else                        % if yes
        hal=holo_arc(end);
    end
    
    if cmplt==1                 % check completeness of the trace
        xholo(:,holostep+2:end)=[];
        aal(holostep+2:end)=[];
    end
    Y=[Y [xholo(:,2:end);aal(2:end)']];
    NYh=NYh+size(aal,1)-1;
    
    if cmplt==1
        break;
    end
    %% numerical continuation for passing through singularity    
    [Z1,Y1,Fail,sgn,cmplt1]=...
        branch_trace_hybrid_4_no_mex(xstart,x,M0,MS,km0,bbt,bb,sgn*svec,...
        param,sgn,aalstart,cont_display,Fail,Y_start,purpose,Tinv,escape,solutionnumber,equationnumber);
    
    Z=[Z Z1];
    Y1(end,:)=Y1(end,:)+aalstart;
    Y = [Y Y1];

    if cmplt1==1
        break;
    elseif cmplt1==0 && Fail(param.nsolu,param.neqt)~=0 && Fail(param.nsolu,param.neqt)~=-1 
        break;
    end
end

[~,Kh]=size(Z);
rpt=0;
Dt=[];
for i=1:Kh-1
     for j=i+1:Kh
         if norm(Z(:,j))~=0 && ( norm(Z(:,j)-Z(:,i),inf)<param.ex || norm(Z(:,j)+Z(:,i),inf)<param.ex)
             Z(:,j)=0;
             rpt=rpt+1;
             Dt(1,rpt)=j;
         end
     end
end
Z(:,Dt)=[];
