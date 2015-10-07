function [AP_tab]=APReadout(AP,a,data_red,kin)

% Creates a cell array containing the affinity parameters for each
% kinase-target-link.

[x,y]=find(a==1);
sc=length(x);
AP_tab=cell(sc,3);
for i=1:sc
    AP_tab{i,1}=data_red(x(i),1);
    AP_tab{i,2}=kin(y(i));
    AP_tab{i,3}=AP(i);
end
end
