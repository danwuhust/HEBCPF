%% Generate quadratic matrices for power flow
clear
clear MakeJacobianD quad_form_vals

%tol=10^(-6);% tolerance
% global bus_n branch_n gen_n numofvar numofcons fac;

%% choose case
fprintf('Obtain the basic data.... \n');
tic;

% mpc=case_illinois200;
% mpc=case118;
% mpc=case57; % 1322 solutions
% mpc=case57mod; % 606 solutions
% mpc=case39; % 176 solutions % 5.3 hrs % 1-core 809 s % parallel 12-core 257.6 s
% mpc=case33bw; % 16 solutions % 160 sec % ellipse basic scaling fa=5, param.mismatch_error=0.1*10^-3
% mpc=case_ieee30; % 472 solutions % 15358.4 s % 1-core 875 s % parallel 12-core 446 s
% mpc=case14; % Numerically unstable. Suggest case14mod
% mpc=case14mod; % 30 solutions % 96.097 s % 1-core 14.86 s % parallel 12-core 9.9 s
% mpc=case14mod2; % reversal load P injection
% mpc=case9Q; % 8 solutions % 14.48 s
% mpc=case9; % 8 solutions % 14.42 s
% mpc=case7Salam; % 4 solutions % 7.9 s
% mpc=case6ww; % 6 solutions % 4.56 s
% mpc=case5loop; % 10 solutions % 5.172 s % ellipse basic scaling fa=5
% mpc=case5Salam; % 10 solutions % 5.078 s % ellipse basic scaling fa=5
% mpc=case5Salam_mod3; % 
% mpc=case4BB0; % 14 solutions % 4.84 s
mpc=case4BBc; % 12 solutions % 3.78 s
% mpc=case4gs; % 6 solutions % 2.55 s
% mpc=case3; % 6 solutions % param.aal_min=2*10^(-5)
% mpc=case3TS; % 6 solutions % 

% increase power injection by a factor
factor = 1+0.19*(0); % maximum = 4/4.2 for 14bus; 2.9/3.5 for 30bus; 2 for 39bus; 1.89/2.35 for 57bus; 1.65/2.0 for 57busmod
mpc.bus(:,3:4)=mpc.bus(:,3:4)*factor;
mpc.gen(:,2:3)=mpc.gen(:,2:3)*factor;
% mpc.bus(:,3)=mpc.bus(:,3)*factor;
% mpc.gen(:,2)=mpc.gen(:,2)*factor;
if factor>=1
    mpc.gen(:,9)=mpc.gen(:,9)*factor;
    mpc.gen(:,4)=mpc.gen(:,9)*0.6;
    mpc.gen(:,5)=-mpc.gen(:,9)*0.6;
end

wait=0;% pause time
[Ybus, ~, ~] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch); % make Y matrix
[bus_n,~]=size(mpc.bus);% bus numbers
numofpv=0;
numofpq=0;
for i=1:bus_n
    if mpc.bus(i,2)==2
        numofpv=numofpv+1;
    end
    if mpc.bus(i,2)==1
        numofpv=numofpq+1;
    end
end
[branch_n,~]=size(mpc.branch);% branch numbers
[gen_n,~]=size(mpc.gen);% generator numbers
numofvar=2*bus_n-1;% number of variables
numofcons=2*bus_n-1;% number of constranints
maxtrial = 1e4;
toc;
%% solution by matpower
tic;
fprintf('\nobtain solution...\n');
result=runpf(mpc);
toc;
%
if result.success == 1
%% obtain quadratic matrices of power balancing equations
tic;
fprintf('\nobtain quadratic matrices of power balancing equations...\n');
[Mp,Mq]=get_quadr_mtrx(Ybus,bus_n);% compute quadratic matrices for power balance equations
% Mps=cell2mat(Mp);
% [MPS(:,1),MPS(:,2),MPS(:,3)]=find(Mps);
% Mqs=cell2mat(Mq);
% [MQS(:,1),MQS(:,2),MQS(:,3)]=find(Mqs);
toc;
%

%% obtain quadratic equations
%
tic;
fprintf('\nobtain quadratic matrices for PV, PQ and slack buses...\n');
[Ma,ba,npg,npd,nq,nv,I]=quadr_matrix(Mp,Mq,mpc,bus_n,gen_n,numofvar);% compute all the quadratic matrices
[solu,solu0]= Solu(result,I.slack,bus_n);% compute an initial solution
Err=zeros(numofcons,1);
for i=1:1:numofcons
   Err(i,1)=solu'*Ma{i}*solu-ba(i); % Err should be around zero for every entry
end
norm(Err,inf) % indicates if everything is calculated correctly
toc;
%}
%% generate ellipses
%
tic;
fprintf('\ngenerate high dimensional ellipses...\n');
[MS,M,bb,M0,km0,bm0,succ,Mss,T] = MakeEllipse(Ma,ba,maxtrial,bus_n,gen_n,numofvar,numofcons); % make ellipses
if succ == 1
% Mss=cell2mat(M);
toc;
%}

%% Resolve initial point
%
tic;
fprintf('\nResolve initial point...\n');
[x,flag] = Resolve_no_mex(MS,bb,M0,km0,solu,numofcons,1); % resolve solution
toc;
%}
%% Tune scaling factor for Jacobian
%
tic;
fprintf('\nTune scaling factor for Jacobian...\n');
[fac] = Tune_Fac(MS,M0,km0,x,numofcons); % tune factor for scaling Jacobian
toc;
%}
else
    fprintf('Fail to generate a base ellipse. Terminated!\n\n');
end

else
    fprintf('Fail to obtain a starting point. Terminated!\n\n');
end

%% generate holomorphic parameters
fprintf('\nGenerate holomorphic parameters...\n');

Tinv=T\eye(numofvar);
degree.max=41;
degree.act=15;
degree.denominator=fix(degree.act/2);
degree.numerator=degree.act-degree.denominator;
