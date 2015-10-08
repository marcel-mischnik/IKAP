function [AP,AP_tab,K,cost,mincost]=FitActivities(a,data_red,kin,ncon,iter)

% Fits the activities of all kinases in kin to the measured data in data_red using 
% the cost functions ComputeCostGlobal.m and ComputeCostLocal.m as well as the truth table a. 

% Ncon corresponds to the number of conditions tested (number of columns 
% with measured values (6 in case of the HeLa data, 8 in case of the insulin data)).
% Iter determines the number of iterations (for reliable results this should be at least 5).

% The output AP is a list of affinity parameters, AP_tab allocates these values to their 
% respective kinase-target-pairs using the function sreadout. K contains the fitted kinase activities 
% as a matrix with ncon columns and length(kin) rows. Cost is a matrix containing all calculated costs, 
% mincost a scalar representing the best optimum found.
% Parameter bounds and starting values should be modified appropriately.

m=size(a,1);
n=size(a,2);
[x,~]=find(a==1);
sc=length(x);
tsm=zeros(sc,iter);
tkm=zeros(n,ncon,iter);
cost=[];
lb1=0.1*ones(sc,1);         
ub1=10*ones(sc,1);
lbk=-15;                % Lower bound for kinase activities. Modify according to your needs.
ubk=15;                 % Upper bound for kinase activities, Modify according to your needs.
lb2=lbk*ones(n,1);   
ub2=ubk*ones(n,1);     
options = optimoptions('fmincon', 'GradObj', 'on');
for i=1:iter
    J1=10000;
    J2=1000;
    ts=10*rand(sc,1);
    tk=15+(-15*rand(n,ncon));   % Starting values for kinase activities. Modify acording to your needs.
    while J1-J2>10
          [ts,J1]=fmincon(@(ts)(ComputeCostGlobal(ts,tk,a,data_red,ncon)), ts, [], [], [], [], lb1, ub1, [], options);
          J=zeros(1,ncon);
          tsm(:,i)=ts;
          for j=3:ncon+2
              tki=tk(:,j-2);
              [tki,J(1,j-2)]=fmincon(@(tki)(ComputeCostLocal(ts,tki,a,data_red(:,j))), tki, [], [], [], [], lb2, ub2, [], options);
              tk(:,j-2)=tki;
              J2=sum(J);
          end
          tkm(:,:,i)=tk;
          cost=[cost; [J1 J2]];
    end
end
J=zeros(1,iter);
for i=1:iter
    [J(i)]=ComputeCostGlobal(tsm(:,i),tkm(:,:,i),a,data_red,ncon);
end
[mincost,g]=min(J);
AP=tsm(:,g);
K=tkm(:,:,g);
[AP_tab]=APReadout(AP,a,data_red,kin);
end
