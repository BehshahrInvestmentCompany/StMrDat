function[data_trans]=DeTransformation(Target, Base_value,Tr)
%[0 1 2 3 4 5]=[no change, Ln, Diff, Double Diff, Diff_Ln, Double Diff_Ln]

[~,N1]=size(Target);
%[~,N2]=size(Exp_Var);

data_adj=[Target(:)];

[T,N]=size(data_adj);
data_trans=NaN(T,N);

for i=1:N
    if Tr(1,i)==0
        data_trans(1:T,i)=data_adj(1:T,i);
    elseif Tr(1,i)==1
        data_trans(1:T,i)=exp(data_adj(1:T,i));
    elseif Tr(1,i)==2
        data_trans(1:T,i)=cumsum(data_adj(1:T,i));
        data_trans(:,i)=data_trans+Base_value*ones(size(data_trans));
    elseif Tr(1,i)==3
        data_trans(1:T,i)=cumsum(cumsum(data_adj(1:T,i)));
        data_trans=data_adj+Base_value*ones(size(data_trans));
    elseif Tr(1,i)==4
        Base_value=log(Base_value);
        data_trans(:,i)=cumsum(data_adj(1:T,i));
        data_trans=data_adj+Base_value*ones(size(data_trans));
        data_trans(:,i)=exp(data_trans(:,i));
    elseif Tr(1,i)==5
        Base_value=log(Base_value);
        data_trans(:,i)=cumsum(cumsum(data_adj(1:T,i)));
        data_trans=data_adj+Base_value*ones(size(data_trans));
        data_trans(:,i)=exp(data_trans(:,i));
    elseif Tr(1,i)==6
        %data_trans(5:T,i)=log(data_adj(5:T,i))-log(data_adj(1:T-4,i));
        
    end
end
% s_data=zeros(1,N);
% for j=1:N
%     if Tr(1,j)==3 || Tr(1,j)==5
%         s_data(1,j)=2;
%     elseif Tr(1,j)==2 || Tr(1,j)==4
%         s_data(1,j)=1;
%     else
%         s_data(1,j)=0;
%     end
% end
% % data_trans=data_trans(max(max(s_data),max(s_data))+1:end,:);
% data_trans=data_trans(max(s_data)+1:end,:);
% Target=data_trans(:,1:N1);
% 
% Exp_Var=data_trans(:,N1+1:N1+N2);
% Dummy=Dummy(max(s_data)+1:end,:);
% Date=Date(max(s_data)+1:end,:);


