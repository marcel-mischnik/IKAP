function [a]=CreateTT(data_red,kcolumn,kin)

% Creates a truth table of logicals (a), indicating if a phosphosphosite is
% phosphorylated by a certain kinase. kcolumn is again the column containing the kinases.

a=zeros(size(data_red,1),size(kin,1));
for i=1:length(data_red)
    for j=1:length(data_red{i,kcolumn})
        kd=data_red{i,kcolumn}{j};
        for l=1:length(kin)
            kk=kin{l};
            f=strcmp(kd,kk);
            if f>0
               a(i,l)=1;
            end
        end
    end
end
end
