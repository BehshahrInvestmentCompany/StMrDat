function [yforecast]=ARXs(Endo, other, dum,one_step,horizon)
%horizon %Forecast horizons
MaxLag=5;
IC=2;    %Information ceriterion 1=AIC, 2=SIC, 3=HQC

% ARDL Model for forecasting 1-4step ahead using Direct method %
% data: inflation  and other predictor variables %
%--------------------------------------------------------------------------
if sum(sum(dum(1:end-horizon,:)))==0
    d=0;
else
    d=1;
    s=sum(dum(1:end-horizon,:),1);
    m=1;
    for i=1:size(s,2)
        if s(1,i)~=0
            cn(1,m)=i;
            m=m+1;
        end
        
    end
end
%-------------------Preallocating Variables--------------------------------
% aic=zeros(5,1);
% sic=zeros(5,1);
% hqc=zeros(5,1);
LagIC=nan(MaxLag,1);
yforecast=zeros(1,horizon);
%--------------------------------------------------------------------------



laginf=lagmatrix(Endo,1:5);
b=cell(5,1);
for i=1:5;
    y=Endo(i+1:end);
    switch d
        case 0
            xx=[ones(length(y),1), laginf(i+1:end, 1:i), other(i+1:end, :)];
        case 1
            xx=[ones(length(y),1), laginf(i+1:end, 1:i), other(i+1:end, :), dum(i+1:end-horizon,cn)];
    end
    b{i,1}=(xx'*xx)\xx'*y;
    yhat=xx*b{i,1};
    e=yhat-y;
    sigma2=(e'*e)/(length(y)-i-size(other,2)-1);
    switch IC
        case 1 %AIC
            LagIC(i,1)=-2*log(sigma2)+2*(i+size(other,2));
        case 2 %SIC
            LagIC(i,1)=log(sigma2)+(i+size(other,2))*log(length(y))/length(y);
        case 3 % HQC
            LagIC(i,1)=log(sigma2)+2*(i+size(other,2))*log(log(length(y)))/length(y);
    end
    
    clear('y', 'yhat', 'e', 'xx') ;
    
end

[ra,~]=find(LagIC==min(LagIC));
% [ra,ca]=find(aic==min(aic));
% [rh,ch]=find(hqc==min(hqc));

% Forecast one-four step ahead %
for h=1:horizon
    laginff=lagmatrix(Endo,[0 1 2 3 4]);
    if d==0
        yforecast(1,h)=[1, laginff(end,1:ra),other(end,:)]*b{ra,1};
    else
        yforecast(1,h)=[1, laginff(end,1:ra),other(end,:), dum(end+h-horizon,cn)]*b{ra,1};
    end
end
if h==1 && isnan(one_step)~=1
    yforecast(1,h)=one_step;
end
Endo=[Endo;yforecast(1,h)];


end


