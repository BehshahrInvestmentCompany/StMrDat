onestep=NaN;
if isnan(onestep)==1
    OS=0;
else
    OS=1;
end
%--- OS=0 if the last quarter is compelete
OneStep=NaN;
%-------------------------Sample Period------------------------------------
fyds=69;    %--- First Year of Data Set
fqds=1;     %--- First quarter of Data Set
%lyds=92;   %--- Last Year of Data Set
%lqds=3;    %--- Last Quarter of Date Set
%-------------------------Forecasting Period for Evaluating Models---------
fyf=89;     %--- First Year of Forecasting
fqf=3;      %--- First Quarter of Forecast
lyf=93;     %--- Last Year of Forecast
lqf=4;      %--- Last Quarter of Forecast
%---------------------Read Data from Excell file---------------------------
Target=xlsread('Data_level.xlsx', 'Target');
Exp_Var=xlsread('Data_level.xlsx', 'Exp_Var');
% F1=xlsread('Data_level.xlsx', 'F1');
% F2=xlsread('Data_level.xlsx', 'F2');
Dummy=xlsread('Data_level.xlsx', 'Dummy');
%---------------------Preaparing Data--------------------------------------
T=[Target(1,2:end), Exp_Var(1,2:end)];
%T is transformation row vector for data
%so that [0 1 2 3 4 5]=[no change, Ln, Diff, Double Diff, Diff_Ln, Double Diff_Ln]
Target=Target(2:end,2:end);
Exp_Var=Exp_Var(2:end,2:end);
% F1=F1(2:end,2:end);
% F2=F2(2:end,2:end);
Date=Dummy(2:end,1);
Dummy=Dummy(2:end,2:end);

[Target Exp_Var Dummy]=Real_time(Target, Exp_Var, Dummy);

f=fullfile(cd,'\IRIS_Tbx');
addpath(f)
irisstartup
[Target_Level_X12, Exp_Var]=deseasonal(Target, Exp_Var);

[Target, Exp_Var, Dummy]=transformation(Target_Level_X12, Exp_Var, Dummy, T);

%--------------------------------------------------------------------------
dl=(length(Date)-4)-length(Target);  %--- Real Time and transformation
sd=4*(fyf-fyds)+(fqf-fqds)+1-dl;     %--- Start Date to Forecasting
nd=4*(lyf-fyf)+(lqf-fqf)+1;          %--- End Date to Forecasting
%--------------------------------------------------------------------------

for i=1:nd+4
    tic
    if OS==1 && i==nd+4
        OneStep=onestep;
    end
    %--------------------------Univariate Models---------------------------
    dum=Dummy(1:sd-4-1+i+4,:);
    Target_=Target(1:sd-4-1+i,1);
   model(i,:,1)=ARdirect(Target_,dum,OneStep);
%     model(i,:,2)=TAR(Target_,OneStep);
%     model(i,:,3)=UM(Target_,OneStep);
%     model(i,:,4)=PureRW(Target_,OneStep);
%     model(i,:,5)=RWDrift(Target_,OneStep);
%     model(i,:,6)=RWAO(Target_,OneStep);
    %-----------------------Multivariate Models----------------------------
    %-------------------------------ARDL-----------------------------------
    Exp_Var_=Exp_Var(1:sd-4-1+i,:);

%     [~,~, ~, ~, f_F1]=factoran(F1(1:sd-4-1+i,:), 3);
%     [~,~, ~, ~, f_F2]=factoran(F2(1:sd-4-1+i,:), 2);
    
    Exp_Var_=[Exp_Var_];
    [~, CC]=size(Exp_Var_);
    
    for j=1:CC
    model(i,:,1+j)=ARDL(Target_,Exp_Var_(:,CC),dum,OneStep);
    end
    %-------------------------------VAR------------------------------------
   % dataset=[Target_, Exp_Var_];
    
 %   [R,C]=size(dataset);
 %   nv=3;                 %%% number of variable in the model %%%
 %   v=nchoosek(2:C,nv);
    
 %   for k=1:size(v,1);
  %      Yraw=[dataset(:,1), dataset(:,v(k,:))];
        %-----------------------------------VAR----------------------------
 %       model(i,:,(1+CC)+k)=BVAR(Yraw,dum,1,OneStep);
        %----------------------------------TVP-VAR-------------------------
%         model(i,:,(6+CC)+length(v)+k)=TVPVAR(Yraw,1,OneStep);
        %---------------------------------ARXs-----------------------------
%         model(i,:,(6+CC)+2*length(v)+k)=ARXs(dataset(:,1),dataset(:,v(k,:)),dum,OneStep);
        %---------------------------------ARXd-----------------------------
        %model(i,:,11+3*length(v)+k)=ARXd(dataset(:,1),dataset(:,v(k,:)),dum,4);
        %---------------------------------ARMAXs---------------------------
        %model(i,:,11+4*length(v)+k)=ARMAXs(dataset(:,1),dataset(:,v(k,:)),dum,4,4);
    end
    %----------------------------------------------------------------------
    
    disp(i)
    tim=toc;
    round(tim*(nd+4-i));
    ['Time remaining to compelete is:' num2str(round(tim*(nd+4-i))) 'Second(s)']

number_models=size(model,3);
%**************************************************************************
%**************************************************************************
psedu_outofsample=zeros(nd,4,number_models);
for i=1:4
psedu_outofsample(:,i,:)=model(5-i:end-i,i,:);
end

outofsample=zeros(4,number_models);
for i=1:4
    outofsample(i,:)=model(end,i,:);
end

%----------------Calculate RMSFE for all Models----------------------------
actual=Target(sd:sd+nd-1,1);
RMSFE1=zeros(4,number_models);
for i=1:4
    for j=1:number_models
        RMSFE1(i,j)=RMSFE(actual, psedu_outofsample(:,i,j));
    end
end


for i=1:4
  [min_RMSFE(i,1), bestmodel(i,1)]=min(RMSFE1(i,:));
end


for i=1:4
    forecast_bestmodel(i,1)=outofsample(i,bestmodel(i,1));
end


%---------------------Foracast Combination(Simple average)-----------------
for i=1:4
    l(i,1)=1;
    for k=1:number_models
        if RMSFE1(i,k)<=RMSFE1(i,1)
            psedu_outofsample_comb(:,i,l(i,1))=psedu_outofsample(:,i,k);
            outofsample_comb(i,l(i,1))=outofsample(i,k);
            l(i,1)=l(i,1)+1;
        end
    end
end

for i=1:4
    psedu_outofsample_combination(:,i)=mean(psedu_outofsample_comb(:,i,1:l(i,1)-1),3);
    outofsample_combination(i,1)=mean(outofsample_comb(i,1:l(i,1)-1),2);
end


for i=1:4
    RMSFE_comb(i,1)=RMSFE(actual, psedu_outofsample_combination(:,i));
end

%----------------Writting forecast of best and combined model to plot fanchart---------

% xlswrite('Result.xlsx', 100*forecast_bestmodel, 'Bestmodel', 'B12:B15');
% xlswrite('Result.xlsx', 100*min_RMSFE, 'Bestmodel', 'D12:D15');
% 
% xlswrite('Result.xlsx', 100*outofsample_combination, 'Combination', 'B7:B10');
% xlswrite('Result.xlsx', 100*RMSFE_comb, 'Combination', 'D7:D10');
% 
% xlswrite('Result.xlsx', Date(end-13:end), 'Calculation', 'a2:a15');
% xlswrite('Result.xlsx', 100*Target(end-9:end), 'Calculation', 'b2:b11');
% xlswrite('Result.xlsx', Target_Level_X12(end-9:end), 'Calculation', 'c2:c11');

save

%construct 2*2 matric from model matrice
for i=1:22
    for j=2:7
        for h=1:4
        Z(i,h)=mean(model(i,:,j));
    end
    end
end


    
    





