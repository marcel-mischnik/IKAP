function [J]=FitIdent(a,pnew,K,lb,ub,data_red,options,sc,c,ncon) 

% Fitting function for identifiability analysis. Should only be used as
% part of IdentKin.m.

J1=100000;
J2=10000;
ki=K;
ki(:,c)=pnew(sc+1:end);
ts=pnew(1:sc);
while J1-J2>10
      [ts,J1]=fmincon(@(ts)(ComputeCostGlobal(ts,ki,a,data_red,ncon)),ts, [], [], [], [], lb(1:sc), ub(1:sc), [], options);
      J=zeros(1,size(K,2));
      for j=3:size(data_red,2)-1
          kp=ki(:,j-2);
          [kin,J(1,j-2)]=fmincon(@(kp)(ComputeCostLocal(ts,kp,a,data_red(:,j))), kp, [], [], [], [], lb(sc+1:end), ub(sc+1:end), [], options);
          ki(:,j-2)=kin;
          J2=sum(J);
      end
end
J=J2;
end
