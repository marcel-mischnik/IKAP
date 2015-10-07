function [psig,pnsig,pDB]=validDB(psig,pnsig,DB)

% Performs a database search with the kinase-target-links in psig and pnsig
% and calculates a p-value based on a fisher exact test. The database
% variable should be a cell array with two columns, the first for targets,
% the second for kinases. Make sure the identifiers conform to those in your data (type, case). 

a=psig(:,1);
b=psig(:,2);
c=DB(:,1);
d=DB(:,2);
e=pnsig(:,1);
f=pnsig(:,2);
for i=1:length(psig)
    st=strcmp(a{i},c);
    sk=strcmp(b{i},d);
    snt=strcmp(e{i},c);
    snk=strcmp(f{i},d);
    for j=1:length(st)
        if st(j)==1 && sk(j)==1
           psig{i,5}=1;
        end
        if snt(j)==1 && snk(j)==1
           pnsig{i,5}=1;
        end
    end
end
su=sum(cell2mat(psig(:,5)));
su1=sum(cell2mat(pnsig(:,5)));
pDB=fexact(su,2*length(psig),su+su1,length(psig));
end
