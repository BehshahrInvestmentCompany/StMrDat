function [yforecast]=ARDL2(Endo_var, Exo_var, dum,one_step, horizon)

% ARDL Model for forecasting 1-4step ahead using Direct method %
% data: inflation  and other predictor variables %
%--------------------------------------------------------------------------
Inx=1:length(Endo_var);
if sum(sum(dum(Inx,:)))==0
    d=0;
else
    d=1;
    s=sum(dum(Inx,:),1);
    m=1;
    for i=1:size(s,2)
        if s(1,i)~=0
            cn(1,m)=i;
            m=m+1;
        end
    end
end
[~,K]=size(Exo_var);
%-------------------Preallocating Variables--------------------------------
% aic=zeros(K,K);
% sic=zeros(K,K);
% hqc=zeros(K,K);
LagIn=nan(K,K);% the index store the best lag by AIC SIC HQC, which is choosen below

yforecast=zeros(1,horizon);
%--------------------------------------------------------------------------
IC=2;    %Information ceriterion 1=AIC, 2=SIC, 3=HQC

for h=1:horizon;
    
    laginf=lagmatrix(Endo_var,h:h+horizon);
    lagother=lagmatrix(Exo_var, h:h+horizon);
    
    b=cell(K,K);
    for i=1:K;
        for j=1:K;
            m=max(h+i,h+j);
            y=Endo_var(m:end);
            switch d
                case 0
                    xx=[ones(length(y),1), laginf(m:end, 1:i), lagother(m:end, 1:j)];
                case 1
                    xx=[ones(length(y),1), laginf(m:end, 1:i), lagother(m:end, 1:j), dum(Inx(m:end),cn)];
            end
            b{i,j}=(xx'*xx)\xx'*y;
            yhat=xx*b{i,j};
            e=yhat-y;
            sigma2=(e'*e)/(length(y)-i-j-1);
            switch IC
                case 1%AIC
                    LagIn(i,j)=-2*log(sigma2)+2*(i+j);
                case 2%SIC
                    LagIn(i,j)=log(sigma2)+(i+j)*log(length(y))/length(y);
                case 3%HQC
                    LagIn(i,j)=log(sigma2)+2*(i+j)*log(log(length(y)))/length(y);
            end
            
            clear('y', 'yhat', 'e', 'xx') ;
            
        end
    end
    %     [rs,cs]=find(sic==min(min(sic)));
    %     [ra,ca]=find(aic==min(min(aic)));
    %     [rh,ch]=find(hqc==min(min(hqc)));
    [ra,ca]=find(LagIn==min(min(LagIn)));
    
    % Forecast one-four step ahead %
    
    laginff=lagmatrix(Endo_var,[0 1 2 3 4 5 6 7 8]);
    lagotherf=lagmatrix(Exo_var, [0 1 2 3 4 5 6 7 8]);
    
    %     switch IC
    %         case 1
    if d==0
        yforecast(1,h)=[1, laginff(end,1:ra),lagotherf(end,1:ca)]*b{ra,ca};
    else
        yforecast(1,h)=[1, laginff(end,1:ra),lagotherf(end,1:ca), dum(end+h-horizon,cn)]*b{ra,ca};
    end
    %         case 2
    %             if d==0
    %                 yforecast(1,h)=[1, laginff(end,1:rs),lagotherf(end,1:cs)]*b{rs,cs};
    %             else
    %                 yforecast(1,h)=[1, laginff(end,1:rs),lagotherf(end,1:cs), dum(end+h-horizon,cn)]*b{rs,cs};
    %             end
    %         case 3
    %             if d==0
    %                 yforecast(1,h)=[1, laginff(end,1:rh),lagotherf(end,1:ch)]*b{rh,ch};
    %             end
    %             yforecast(1,h)=[1, laginff(end,1:rh),lagotherf(end,1:ch), dum(end+h-horizon,cn)]*b{rh,ch};
end
if isnan(one_step)~=1
    yforecast(1,1)=one_step;
end
end





