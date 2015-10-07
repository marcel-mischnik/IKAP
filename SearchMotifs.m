function [data]=SearchMotifs(data,PSP)

% Searches the PhosphoSitePlus database (download Kinase-Substrate-Dataset.gz 
% at http://www.phosphosite.org/staticDownloads.do) for kinases that are 
% known to phosphorylate the central amino acid of the given motif. The motifs 
% can have an arbitrary length, however all motifs should be equally long. 
% The function adds a further column to the data containing the respective kinase(s).

% PSP should be a two-column cell array generated from the downloaded file with 
% kinase names (e.g. gene names) in the first column and motif in the second (no headers). 
% The motifs in PSP are assumed to be 15 amino acids long.

% Data should be a cell array beginning with one column for the protein name (or ID) 
% the phosphosite is located on, followed by one column for the motif and a 
% variable number of numeric columns. All IKAP functions assume a total number 
% of two annotation columns. If more annotation columns shall be included, the functions 
% have to be altered where neccessary. The header row should be removed.

data=[data cell(size(data,1),1)];
le=length(data{2,2});
mdata=(le+1)/2;
mpsp=8;
if isnumeric(data{1,3})
    start=1;
else
    start=2;
end
if mdata>=mpsp
    for i=start:length(data)
        motdata=data{i,2}(mdata-7:mdata+7);
        s=strcmp(motdata,PSP(:,2));
        if sum(s)>0
            x=find(s==1);
            for j=1:length(x)
                if isempty(data{i,end}) && ~isempty(PSP{x(j),1})
                   data{i,end}=PSP{x(j),1};
                elseif ~isempty(data{i,end}) && ~isempty(PSP{x(j),1})
                   data{i,end}=[data{i,end} PSP(x(j),1)];
                end
            end
        end
    end
else
    for i=1:length(PSP)
        motpsp=PSP{i,2}(8-mdata:8+mdata);
        s=strcmp(motpsp,data(:,2));
        if sum(s)>0
            x=find(s==1);
            for j=1:length(x)
                if isempty(data{x(j),end})
                   data{x(j),end}=PSP{i,1};
                else 
                   data{x(j),end}=[data{x(j),end} PSP{i,1}];
                end
            end
        end
    end
end

for h=start:length(data)
    data{h,end}=upper(data{h,end});
    if ischar(data{h,end})
        data{h,end}=cellstr(data{h,end});
    end
    data{h,end}=unique(data{h,end});
end
end
