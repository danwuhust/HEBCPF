%runbook.  Script to trace and save traces with bookkeeping.
%Initialize:  start with solution x.

% global bus_n gen_n numofcons numofvar;
Tholo=0;
Tresol=0;
Tcont=0;
Tdata=0;
totalTime=0;
if exist('NY','var')
    
else
    NY=0;
    NYh=0;
end
clear Ysave;

tcpu_start = cputime; % recording the starting time
numberofequations = length(bb);
% solutionnumber = 1;
equationnumber=1;
N.PV = find(mpc.bus(:,2)==2);
N.PQ = find(mpc.bus(:,2)==1);
N.slack = find(mpc.bus(:,2)==3);
%{
if size(N.PQ,1)>=2
    if N.PQ(1)<N.slack && N.PQ(2)<N.slack
        variablelist = [N.PQ(1), N.PQ(2)] + bus_n; 
    elseif N.PQ(1)<N.slack && N.PQ(2)>N.slack
        variablelist = [N.PQ(1)+bus_n, N.PQ(2) + bus_n-1];
    elseif N.PQ(1)>N.slack && N.PQ(2)<N.slack
        variablelist = [N.PQ(1)+bus_n-1, N.PQ(2) + bus_n];
    elseif N.PQ(1)>N.slack && N.PQ(2)>N.slack
        variablelist = [N.PQ(1), N.PQ(2)] + bus_n -1; 
    end
else
    variablelist = bus_n+1:2*bus_n-1;
end
%}
variablelist = 1:numofvar;
numberofvariables = length(variablelist);
Fail=[]; % failure recording book
count = 1; % count the number of traces

Trace=cell(1,1); % trace record
initbypass =0; % if it is 0, starting from the very beginning; if it is 1, starting from the last solution

%% if VBook exists, check if VBook matches Zsave, and start from the existing solution
if exist('VBook','var')
    numberofsolutions = size(VBook,1)
    if size(Zsave,2)==numberofsolutions % Vbook record matches Zsave record
        initbypass = 1; % starting from the last solution
        count = max(max(VBook))+1 % count the latest number of traces and proceed to the next one
    end
end

%% 
if initbypass==1
%     solutionnumber = 1;
else % initialize, starting from the first solution and the first trace
    clear VBook Zsave Vsave;
    
    Zsave=[];
    Zsave_chunk = 100;  % chunk preallocation size
    Vsave=[];
    
    svec = zeros(size(bb)); 
    svec(1) = 1;
        
    solutionnumber = 1;

    tic;
%     traceloops_20181012; % first path
%     hybrid_traceloops_2;
%     hybrid_traceloops_3;
    hybrid_traceloops_4_with_mex;

    Time=toc;
    NY=NY+size(Y,2);
    totalTime=totalTime+Time;
%     pause(3)
%     Trace{1,1}=Y;
    
    numberofsolutions = size(Z,2);
        
    VBook(1:numberofsolutions,1) = count;  %bookkeeping. 
    
    n_rows = size(Z, 1);
    Zsave = zeros(n_rows, max(numberofsolutions, Zsave_chunk));
    Zsave(:, 1:numberofsolutions) = Z;
%     count = count + 1; 
    fprintf('No. Solu: %d, Solu No.: %d, No. Eqt: %d, Eqt No.: %d, Time: %g \n',numberofsolutions, solutionnumber, numberofequations, 1,Time);
    pause(wait);
end

%%
while numberofsolutions+1>solutionnumber
    equationqueue = setdiff(1:numberofequations,find(VBook(solutionnumber,:)>0));   
    Neqnqueue = size(equationqueue,2); % length of equation queue
    New_solu = cell(Neqnqueue,1); % new solution cell
    New_solu_match = cell(Neqnqueue,1);
    Zsave_ind = 1:numberofsolutions;

    %% trace different loops
    x_ini = Zsave(:,solutionnumber);
    for eqnN = 1:Neqnqueue %equationnumber = equationqueue
        warning off MATLAB:nearlySingularMatrix;
        equationnumber = equationqueue(eqnN);
        x = x_ini;
        Fail = -1;
        svec = zeros(numberofequations,1); 
        svec(equationnumber) = 1;
        
        tic;
%         traceloops_20181012;
%         hybrid_traceloops_2;
%         hybrid_traceloops_3;
        hybrid_traceloops_4_with_mex;

%%
        Time=toc;
        NY=NY+size(Y,2);
        totalTime=totalTime+Time;

        New_solu{eqnN,1}=Z;

        fprintf('No. Solu: %d, Solu No.: %d, No. Eqt: %d, Eqt No.: %d, Time: %g\n',numberofsolutions, solutionnumber, numberofequations, equationnumber,Time);
        pause(wait);
    end
    
    %% collecting new solutions
    for eqind0 = 1:Neqnqueue
        count = count + 1;
        New_solu_match{eqind0}(1,1)=solutionnumber;
        for eqind1 = 2:size(New_solu{eqind0},2)
            dZoldnew = sum(abs(Zsave(variablelist,1:numberofsolutions)-New_solu{eqind0}(variablelist,eqind1)),1);
            eqind2 = find(dZoldnew<4*10^(-7),1);

            dZoldnew2 = sum(abs(Zsave(variablelist,1:numberofsolutions)+New_solu{eqind0}(variablelist,eqind1)),1);
            eqind3 = find(dZoldnew2<4*10^(-7),1);
            if ~isempty(eqind2)
                New_solu_match{eqind0}(1,eqind1)=eqind2;
            elseif ~isempty(eqind3)
                New_solu_match{eqind0}(1,eqind1)=eqind3;
            else
                numberofsolutions = numberofsolutions+1;
                if numberofsolutions > size(Zsave, 2)
                    new_capacity=max(numberofsolutions,2*size(Zsave,2));
                    Zsave(:,new_capacity)=0;
                end
                Zsave(:, numberofsolutions) = New_solu{eqind0}(:, eqind1);
                New_solu_match{eqind0}(1,eqind1)=numberofsolutions;

            end
        end
        VBook(New_solu_match{eqind0},equationqueue(eqind0))=count;
    end

    solutionnumber = solutionnumber +1;
end
Zsave = Zsave(:, 1:numberofsolutions);  % trim the extra preallocated columns
[Zsave,VBook]=canonicalize_solutions(Zsave,VBook);
tcpu_over=cputime;
fprintf('Overall executing time: %g\n',tcpu_over-tcpu_start);
