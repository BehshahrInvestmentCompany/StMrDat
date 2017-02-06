function [Exp_Var]=deseasonal(Exp_Var)
%[T1,N1]=size(Target);
%[T2,N2]=size(Exp_Var);

data=[Exp_Var];
addpath current_path\IRIS_Tbx; irisstartup
%addpath current_path\IRIS_Tbx; irisstartup
q=qq(1900,1):qq(2000,1);
q=q';


[T N]=size(data);
for i=1:N
    x1=tseries(q(1:T,1),data(1:T,i));
    x2=x12(x1);
    x3=x2.data;
    data_Adj(:,i)=x3;
end
%Target=data_Adj(:,1:N1);
Exp_Var=data_Adj(:,:);

end


