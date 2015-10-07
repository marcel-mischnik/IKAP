function [plJ,plk]=IdentKin(ks,a,c,AP,K,data_red,ncon)

% Calculates identifiability profiles for the kinases given in ks. This
% should be a vector containing the kinase numbers (not names).

% Inputs:
% ks: vector with kinase numbers
% a: truth table 
% c: column number of the condition in data_red the identifiabiliy shall be calculated on
% AP: vector with fitted affinity parameters
% K: matrix with fitted kinase activities
% data_red: reduced dataset

% Outputs:
% plJ: matrix with calculated costs along the dimension of the respective kinase
% plk: matrix with tested kinase values

p=[AP; K(:,c-2)];
plJ=zeros(length(ks),7);
plk=zeros(length(ks),7);
adds=-3:3;
[x,~]=find(a==1);
sc=length(x);
lkpr=size(a,2);
lbk=-15;            % Lower bound for kinase activities. Modify according to your needs.
ubk=15;             % Upper bound for kinase activities. Modify according to your needs.
options = optimoptions('fmincon', 'GradObj', 'on');
parfor l=1:length(ks)
    i=ks(l);
    for j=1:7
        pnew=p;
        pnew(sc+i)=p(sc+i)+adds(j);
        plk(l,j)=p(sc+i)+adds(j);
        if j==4
            [plJ(l,j)]=ComputeCostGlobal(AP,K,a,data_red,ncon);
        else
            lb=[0.1*ones(sc,1); lbk*ones(lkpr,1)];
            ub=[10*ones(sc,1); ubk*ones(lkpr,1)];
            lb(sc+i)=pnew(sc+i);
            ub(sc+i)=pnew(sc+i);
            plJ(l,j)=FitIdent(a,pnew,K,lb,ub,data_red,options,sc,c,ncon);
        end
    end    
end
end
