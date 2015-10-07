function [dist,valids]=ComputeDistances(data,ncon,kin,K)

% Calculates p-values for each kinase-phosphosite-pair by means of a
% correlation coefficient.

% Inputs:
% data: complete data matrix
% ncon: number of conditions
% kin: kinase list
% K: estimated kinase activities

% Outputs:
% dist: cell array allocating a p-value to each combination of 
% phosphosites (y-axis) and kinases (x-axis).
% valids: numeric matrix indicating on how many valid measurement values 
% the correlation coefficient is calculated in each case.

dist=cell(size(data_red,1),size(kin,1)+2);
valids=zeros(size(data_red,1),1);

for i=1:size(data_red,1)
    da=data_red(i,1:size(data_red,2));
    da(strcmp('NaN',da))={NaN};
    dist(i,1:2)=da(1:2);
    dat=cell2mat(da(3:2+ncon));
    d=find(~isnan(dat));
    for m=1:length(kin)
        [~,P]=corrcoef(dat(d),K(m,d));
        minp=min(P);
        if ~isnan(minp(1))
            dist{i,m+2}=minp(1);
        else
            dist{i,m+2}=1;
        end
        valids(i)=length(d);
    end
end
dist=[cell(1,length(kin)+2); dist];
dist(1,3:end)=kin';
end
