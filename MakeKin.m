function [data_red,kin]=MakeKin(data,kcolumn,kproteome)

% Produces a list of kinases (kin) that are found by SearchMotifs.m.
% kcolumn corresponds to the column number of the data matrix containing the kinases.
% 
% If proteome data is available (optional one column cell array containing protein IDs (kproteome)), 
% the list is reduced to those kinases which are present in the proteome.  

% Generates a new data array (data_red) that only includes those sites for which a kinase is known. 

kin=cell(1,1);
for i=1:length(data)
    kin=[kin data{i,kcolumn}];
end
kin(1)=[];
kin=unique(kin);

if nargin==3
    kpr2={};
    for i=2:length(kin)
        k1=kin{i};
        for j=1:length(kproteome)
            k2=kproteome{j};
            f=strcmp(k1,k2);
            if sum(f)>0
                kpr2=[kpr2 k1];
            end
        end
    end
    %kpr2(1)=[];
    kpr2=unique(kpr2);
    kin=kpr2;
end
if nargin==2
    kin(1)=[];
end
data_red=cell(1,kcolumn);
for j=2:length(data)
    b=0;
    if ~isempty(data{j,kcolumn})
        for i=1:length(data{j,kcolumn})
            if sum(strcmp(kin(:),data{j,kcolumn}{i}))>0
               b=1;
            end
        end
        if b==1
           data_red=[data_red; data(j,:)];
        end
    end
end
data_red(1,:)=[];
kin=kin';
end
