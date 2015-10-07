function [psig,pnsig,pSUM,pNUM]=validMOT(psig,pnsig,nwk)

% Searches for motifs in the NetworKIN database and calculates p-values
% based on the respective likelihoods using a fisher exact test. The
% database can be downloaded at http://www.networkin.info/donwload.shtml.
% pSUM is computed based on the sum of the likelihoods in psig and pnsig, 
% whereas pNUM refers to the number of likelihoods larger
% than 1. The variable nwk should be a cell array with four columns:
% targets, likelihoods, kinases, motifs.

psig(:,6)={0};
pnsig(:,6)={0};
f=nwk(:,4);
g=psig(:,3);
h=pnsig(:,3);
bs=cell(length(nwk),1);
b1s=cell(length(nwk),1);
le=length(psig{1,3});
les=(le-1)/2;

parfor i=1:length(nwk)
    s1=strfind(g(:),f{i});
    s2=strfind(h(:),f{i});
    [x1]=cellfun(@isempty,s1);
    [bs{i}]=find(x1==0);
    [x2]=cellfun(@isempty,s2);
    [b1s{i}]=find(x2==0);
end
for i=1:length(bs)
    if ~isempty(bs{i})
       for j=1:length(bs{i})
           if strcmp(psig{bs{i}(j),2},nwk{i,1}) && strcmp(psig{bs{i}(j),1},nwk{i,3})
              if strcmp(nwk{i,4}(1:5),psig{bs{i}(j),3}(les-4:les)) || strcmp(nwk{i,4}(7:end),psig{bs{i}(j),3}(les+1:les+5))
                 psig{bs{i}(j),6}=nwk{i,2};
              end
           end
       end
    end
    if ~isempty(b1s)
       for j=1:length(b1s{i})
           if strcmp(pnsig{b1s{i}(j),2},nwk{i,1}) && strcmp(pnsig{b1s{i}(j),1},nwk{i,3})
              if strcmp(nwk{i,4}(1:5),pnsig{b1s{i}(j),3}(les-4:les)) || strcmp(nwk{i,4}(7:end),pnsig{b1s{i}(j),3}(les+1:les+5))
                 pnsig{b1s{i}(j),6}=nwk{i,2};
              end
           end
       end
    end    
end

m1=floor(sum(cell2mat(psig(:,6))));
m2=floor(sum(cell2mat(pnsig(:,6))));
sigmot1=length(find(cell2mat(psig(:,6))>=1));
sigmot2=length(find(cell2mat(pnsig(:,6))>=1));
pSUM=fexact(m1,2*length(psig),m1+m2,length(psig));
pNUM=fexact(sigmot1,2*length(psig),sigmot1+sigmot2,length(psig));
end
