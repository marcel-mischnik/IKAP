function [J,grad]=ComputeCostLocal(ts,tki,a,b)

% Calculates the cost and the kinase activity gradients for the condition
% given in b (one data column of data_red) using the affinity parameters
% ts, the local kinase activities tki and the truth table a.

m=size(a,1);
n=size(a,2);
J=0;
grad=zeros(n,1);
p=[ts; tki];
p1=zeros(m,n);
[x,y]=find(a==1);
sc=length(x);

for i=1:sc
    p1(x(i),y(i))=p(i);
end
p1=[p1; p(sc+1:end)'];

for i=1:m
    if ~isnan(b{i})
       a_ex=p1(m+1,:).*p1(i,:);
       avg=sum(a_ex)/sum(p1(i,:));
       if ~isnan(avg)
           J=J+(avg-b{i})^2;
       end
    end
end

for l=1:n
    N=0;
    for i=1:m
        if ~isnan(b{i}) && a(i,l)~=0
           a_ex=p1(m+1,:).*p1(i,:);
           avg=sum(a_ex)/sum(p1(i,:));
           if ~isnan(avg)
               N=N+(2*(avg-b{i})*(p1(i,l)/sum(p1(i,:))));
           end
        end
    end
    grad(l)=N;
end
end
