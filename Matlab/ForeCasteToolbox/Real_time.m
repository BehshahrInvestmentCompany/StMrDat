function [Target ,Exp_Var, Dummy,Date]=Real_time(Target, Exp_Var,  Dummy,Date)
% it must be revised
% The Goal of this function is to remove the extra date which has no value
% [T1,N1]=size(Exp_Var);
inx=any(isnan([Target, Exp_Var,  Dummy]),2);
% data=[Exp_Var];
% [T,N]=size(data);
% index=zeros(1,N);
% for i=1:N
%     for j=1:T
%         if isnan(data(j,i))==1
%             index(1,i)=index(1,i)+1;
%         else
%             index(1,i)=index(1,i)+0;
%         end
%     end
% end
% 
% for i=1:N
%     data(1+index(1,i):end,i)=data(1:end-index(1,i),i);
% end
% data=data(max(index)+1:end,:);
Exp_Var(inx,:)=[];
Target(inx,:)=[];
Dummy(inx,:)=[];


Date(inx,:)=[];
