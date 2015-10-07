function [J,grad]=ComputeCostGlobal(ts,tk,a,data_red,ncon)

% Calculates the cost and the gradients for the affinity parameters over 
% all conditions using the affinity factors (ts), the kinase activities (tk), 
% the truth table (a) and the data (data_red). 

% Ncon corresponds to the number of conditions measured (number of data
% columns).

m=size(a,1);
n=size(a,2);
J=0;
grad=[];
p1=zeros(m,n);
[x,y]=find(a==1);
for i=1:length(x)
    p1(x(i),y(i))=ts(i);
end
for i=1:ncon
    p1=[p1; tk(:,i)'];
end

for j=1:ncon
    for i=1:m
        if ~isnan(data_red{i,j+2})
           a_ex=p1(m+j,:).*p1(i,:);
           avg=sum(a_ex)/sum(p1(i,:));
           if ~isnan(avg)
               J=J+(avg-data_red{i,j+2})^2;
           end
        end
    end
end

for i=1:m
    for j=1:n
        g=0;
        for l=1:ncon
            if ~isnan(data_red{i,l+2}) && a(i,j)~=0
               g=g+(2*((p1(m+l,j)*p1(i,j))/sum(p1(i,:))-data_red{i,l+2})*(p1(m+l,j)*sum(p1(i,:))-sum(p1(m+l,:).*p1(i,:)))/(sum(p1(i,:)))^2);
            end
        end
        if a(i,j)~=0
           grad=[grad; g];
        end
    end
end

end
