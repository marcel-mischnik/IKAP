function [pnsig]=MakePnsig(dist,psig,kin)

% Produces a list of randomly selected links (pnsig) that has the same size as
% psig.

pnsig=cell(length(psig),5);
r1=randi(size(dist,1)-1,length(psig),1);
r2=randi(length(kin),length(psig),1);
cd=dist(2:end,3:end);
sites=dist(2:end,1:2);
pnsig(:,1)=sites(r1,1);
pnsig(:,2)=kin(r2,1);
pnsig(:,3)=sites(r1,2);
for i=1:length(psig)
    pnsig(i,4)=cd(r1(i),r2(i));
    if isempty(pnsig{i,1}) || isnumeric(pnsig{i,1})
        pnsig{i,1}='NAN';
    end
end
pnsig(:,5)={0};
pnsig(:,1)=upper(pnsig(:,1));
end
