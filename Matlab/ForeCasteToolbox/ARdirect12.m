function [yforecast]=ARdirect(inflation,dum,one_step)
% AR(p) Forecasting inflation bye direct approach %
%--------------------------------------------------------------------------
if sum(sum(dum(1:end-4,:)))==0
    d=0;
else
    d=1;
    s=sum(dum(1:end-4,:),1);
    m=1;
    for i=1:size(s,2)
        if s(1,i)~=0
            cn(1,m)=i;
            m=m+1;
        end
        
    end
end

%--------------------Preallocating Variables-------------------------------
beta=cell(5,1);
aic=zeros(5,1);
sic=zeros(5,1);
hqc=zeros(5,1);
yforecast=zeros(1,4);
%--------------------------------------------------------------------------
IC=2;    %Information ceriterion 1=AIC, 2=SIC, 3=HQC

for h=1:4;
    for i=1:5;
        laginf=lagmatrix(inflation, h:h+i-1);
        y=inflation(h+i:end);
        switch d
            case 0
                xx=[ones(length(y),1), laginf(h+i:end,:)];
            case 1
                xx=[ones(length(y),1), laginf(h+i:end,:), dum(h+i:end-4,cn)];
        end
        beta{i,1}=xx\y;
        yhat=xx*beta{i,1};
        e=y-yhat;
        sigma2=(e'*e)/(length(y)-i-1);
        switch IC
            case 1
                aic(i,1)=-2*log(sigma2)+2*(i);
            case 2
                sic(i,1)=log(sigma2)+(i)*log(length(y))/length(y);
            case 3
                hqc(i,1)=log(sigma2)+2*(i)*log(log(length(y)))/length(y);
        end
        
        clear y xx yhat e
    end
    
    [min_aic ra]=min(aic);
    [min_sic rs]=min(sic);
    [min_hqc rh]=min(hqc);
    
    switch IC
        
        case 1
            laginff=lagmatrix(inflation, 0:ra-1);              % shif data for forecasting one-step ahead %
            if d==0
                yforecast(1,h)=[1, laginff(end,:)]*beta{ra,1};
            else
                yforecast(1,h)=[1, laginff(end,:), dum(end+h-4,cn)]*beta{ra,1};
            end
            
        case 2
            laginff=lagmatrix(inflation, 0:rs-1);              % shif data for forecasting one-step ahead %
            if d==0
                yforecast(1,h)=[1, laginff(end,:)]*beta{rs,1};
            else
                yforecast(1,h)=[1, laginff(end,:), dum(end+h-4,cn)]*beta{rs,1};
            end
            
        case 3
            laginff=lagmatrix(inflation, 0:rh-1);              % shif data for forecasting one-step ahead %
            if d==0
                yforecast(1,h)=[1, laginff(end,:)]*beta{rh,1};
            else
                yforecast(1,h)=[1, laginff(end,:), dum(end+h-4,cn)]*beta{rh,1};
            end
    end
    clear laginff
    if isnan(one_step)~=1
    yforecast(1,1)=one_step;
    end
end

