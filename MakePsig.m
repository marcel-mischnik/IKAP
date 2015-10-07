function [psig]=MakePsig(dist,kin,q)

% Produces a list of significant links applying the Benjamini-Hochberg procedure 
% with a desired false discovery rate q. 
% The resulting cell array contains one column with protein IDs, one with 
% kinase IDs and one with the respective p-value.

[~,crit_p]=fdr_bh(cell2mat(dist(2:end,3:end)),q);

if iscell(dist)
    di=cell2mat(dist(2:end,3:end));
else
    di=dist;
end
[x,y]=find(di<=crit_p);
psig=cell(length(x),5);
psig(:,1)=dist(x,1);
psig(:,2)=kin(y);
psig(:,3)=dist(x,2);
for i=1:length(psig)
    psig{i,4}=di(x(i),y(i));
    if isnumeric(psig{i,1}) || isempty(psig{i,1})
        psig{i,1}='NAN';
    else
        psig{i,1}=upper(psig{i,1});
    end    
end
psig(:,5)={0};
end
