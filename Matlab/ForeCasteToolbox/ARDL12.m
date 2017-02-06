function [yforecast]=ARDL(inflation, other, dum,one_step)

% ARDL Model for forecasting 1-4step ahead using Direct method %
% data: inflation  and other predictor variables %
%--------------------------------------------------------------------------
if sum(sum(dum(1:end-12,:)))==0
    d=0;
else
    d=1;
    s=sum(dum(1:end-12,:),1);
    m=1;
    for i=1:size(s,2)
        if s(1,i)~=0
            cn(1,m)=i;
            m=m+1;
        end
        
    end
end
%-------------------Preallocating Variables--------------------------------
aic=zeros(5,5);
sic=zeros(5,5);
hqc=zeros(5,5);
yforecast=zeros(1,12);
%--------------------------------------------------------------------------
IC=2;    %Information ceriterion 1=AIC, 2=SIC, 3=HQC

for h=1:12;
    
    laginf=lagmatrix(inflation,h:h+12);
    lagother=lagmatrix(other, h:h+12);
    b=cell(5,5);
    for i=1:5;
        for j=1:5;
            m=max(h+i,h+j);
            y=inflation(m:end);
            switch d
                case 0
                    xx=[ones(length(y),1), laginf(m:end, 1:i), lagother(m:end, 1:j)];
                case 1
                    xx=[ones(length(y),1), laginf(m:end, 1:i), lagother(m:end, 1:j), dum(m:end-12,cn)];
            end
            b{i,j}=(xx'*xx)\xx'*y;
            yhat=xx*b{i,j};
            e=yhat-y;
            sigma2=(e'*e)/(length(y)-i-j-1);
            switch IC
                case 1
                    aic(i,j)=-2*log(sigma2)+2*(i+j);
                case 2
                    sic(i,j)=log(sigma2)+(i+j)*log(length(y))/length(y);
                case 3
                    hqc(i,j)=log(sigma2)+2*(i+j)*log(log(length(y)))/length(y);
            end
            
            clear('y', 'yhat', 'e', 'xx') ;
            
        end
    end
    [rs,cs]=find(sic==min(min(sic)));
    [ra,ca]=find(aic==min(min(aic)));
    [rh,ch]=find(hqc==min(min(hqc)));
    
    % Forecast one-four step ahead %
    
    laginff=lagmatrix(inflation,[0 1 2 3 4 5 6 7 8 9 10 11 12]);
    lagotherf=lagmatrix(other, [0 1 2 3 4 5 6 7 8 9 10 11 12]);
    
    switch IC
        case 1
            if d==0
                yforecast(1,h)=[1, laginff(end,1:ra),lagotherf(end,1:ca)]*b{ra,ca};
            else
                yforecast(1,h)=[1, laginff(end,1:ra),lagotherf(end,1:ca), dum(end+h-12,cn)]*b{ra,ca};
            end
        case 2
            if d==0
                yforecast(1,h)=[1, laginff(end,1:rs),lagotherf(end,1:cs)]*b{rs,cs};
            else
                yforecast(1,h)=[1, laginff(end,1:rs),lagotherf(end,1:cs), dum(end+h-12,cn)]*b{rs,cs};
            end
        case 3
            if d==0
                yforecast(1,h)=[1, laginff(end,1:rh),lagotherf(end,1:ch)]*b{rh,ch};
            end
            yforecast(1,h)=[1, laginff(end,1:rh),lagotherf(end,1:ch), dum(end+h-12,cn)]*b{rh,ch};
    end
end

if isnan(one_step)~=1
    yforecast(1,1)=one_step;
end
end


