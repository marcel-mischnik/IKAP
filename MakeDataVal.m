function [data_val]=MakeDataVal(data)

% Creates a dataset from which the phosphosites that were used for
% parameter estimation are removed. This guerantees for an unbiased validation. 

data_val={};
for i=1:length(data)
    if isempty(data{i,end})
        data_val=[data_val; data(i,:)];
    end
end
end
