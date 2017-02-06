function [forecast]=TAR(inf, one_step,horizon)
% It must be Revised
MaxLag=5;
IC=2; %Information ceriterion 1=AIC(defalut), 2=SIC, 3=HQC

sortinf=sort(inf);
cmin=round(0.4*length(sortinf));
cmax=round(0.6*length(sortinf));
cc=sortinf(cmin:cmax,1);

%---------------Preallocating Variables------------------------------------
sizecc=length(cc);
beta1=cell(sizecc,1);
beta2=cell(sizecc,1);
SIC1=zeros(sizecc,1);
SIC2=zeros(sizecc,1);
SIC=zeros(sizecc,1);
%--------------------------------------------------------------------------


for i=1:length(cc)
    laginf=lagmatrix(inf, 0:MaxLag);
    laginf=laginf(MaxLag+1:end,:);
    m=1;
    n=1;
    for j=1:length(laginf)
        if laginf(j,2)<=cc(i)
            y1(m,:)=laginf(j,:);
            m=m+1;
        else
            y2(n,:)=laginf(j,:);
            n=n+1;
        end
    end
    [beta1{i,1}, SIC1(i,1)]=TAR_AR(y1,IC,MaxLag);
    [beta2{i,1}, SIC2(i,1)]=TAR_AR(y2,IC,MaxLag);
    SIC(i,1)=SIC1(i,1)+SIC2(i,1);
    
    clear y1 y2
end

rs=find(SIC==min(min(SIC)));
if length(rs)>=1
    rs=rs(1);
end
%-------------------------OneStep------------------------------------------
laginff=lagmatrix(inf, 0:4);
if laginff(end, 1)<=cc(rs)
    onestep=[1, laginff(end, 1:length(beta1{rs,1})-1)]*beta1{rs,1};
else
    onestep=[1, laginff(end, 1:length(beta2{rs,1})-1)]*beta2{rs,1};
end
clear laginff
if isnan(one_step)~=1
    onestep=one_step;
end
%-------------------------TwoSteps-----------------------------------------
inf=[inf;onestep];
laginff=lagmatrix(inf, 0:4);
if laginff(end, 1)<=cc(rs)
    twostep=[1, laginff(end, 1:length(beta1{rs,1})-1)]*beta1{rs,1};
else
    twostep=[1, laginff(end, 1:length(beta2{rs,1})-1)]*beta2{rs,1};
end
clear laginff
%-------------------------ThreeSteps---------------------------------------
inf=[inf;twostep];
laginff=lagmatrix(inf, 0:4);
if laginff(end, 1)<=cc(rs)
    threestep=[1, laginff(end, 1:length(beta1{rs,1})-1)]*beta1{rs,1};
else
    threestep=[1, laginff(end, 1:length(beta2{rs,1})-1)]*beta2{rs,1};
end
clear laginff
%-------------------------fourSteps----------------------------------------
inf=[inf;threestep];
laginff=lagmatrix(inf, 0:4);
if laginff(end, 1)<=cc(rs)
    fourstep=[1, laginff(end, 1:length(beta1{rs,1})-1)]*beta1{rs,1};
else
    fourstep=[1, laginff(end, 1:length(beta2{rs,1})-1)]*beta2{rs,1};
end
clear laginff

forecast=[onestep, twostep, threestep, fourstep];

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%------------------Implicit Function(s)------------------------------------
function [beta, SIC]=TAR_AR(y,IC,MaxLag)
% y : Data
% IC : Information ceriterion 1=AIC, 2=SIC, 3=HQC
if nargin<2
    IC=1;
end
if nargin<3
    MaxLag=5;
end
%--------------------Preallocating Variables-------------------------------
b=cell(MaxLag,1);
aic=zeros(MaxLag,1);
sic=zeros(MaxLag,1);
hqc=zeros(MaxLag,1);
%--------------------------------------------------------------------------
for j=1:MaxLag;
    
    yy=y(:,1);
    xx=[ones(length(yy),1), y(:,2:1+j)];
    b{j,1}=(xx'*xx)\xx'*yy;
    yhat=xx*b{j,1};
    e=yy-yhat;
    sigma2=(e'*e)/(length(yy)-j-1);
    aic(j,1)=-2*log(sigma2)+2*(j);
    sic(j,1)=log(sigma2)+(j)*log(length(yy))/length(yy);
    hqc(j,1)=log(sigma2)+2*(j)*log(log(length(yy)))/length(yy);
    
    clear laginf yy xx yhat e sigma2
end;
ra=find(aic==min(min(aic)));
rs=find(sic==min(min(sic)));
rh=find(hqc==min(min(hqc)));
switch  IC
    case 2
        beta=b{rs,1};
        SIC=sic(rs,1);
    case 3
        beta=b{ra,1};
        SIC=sic(ra,1);
    otherwise
        beta=b{rh,1};
        SIC=sic(rh,1);
end


