
%{
f=fullfile(cd,'\IRIS_Tbx');
addpath(f)
irisstartup
%}
% Version 2
% Base on Karami and Bayat
% Modified By P.Davoudi: Pedram.Davoudi@gmail.com
%% Initializing
clear
clc
%****************************
onestep=NaN;
horizon=4;
% if isnan(onestep)==1
%     OS=0;
% else
%     OS=1;
% end
%--- OS=0 if the last quarter is compelete
OneStep=nan;% OneStep=NaN;

% %---------------------Read Data from Excell file---------------------------
% Target=xlsread('Input\Data_level.xlsx', 'Target');
% Exp_Var=xlsread('Input\Data_level.xlsx', 'Exp_Var');
% % F1=xlsread('Data_level.xlsx', 'F1');
% % F2=xlsread('Data_level.xlsx', 'F2');
% Dummy=xlsread('Input\Data_level.xlsx', 'Dummy');
Data=dataset('xls','Input\Data.xlsx','sheet','Data');
Target_var_name={'CPI'};
Dum_var_name={'D1','D2'};

T_num=1;
% ---------------------------------------------------
Var_names= Data.Properties.VarNames;
Target_var_Pos=strcmp(Var_names,Target_var_name{T_num});
Target=double(Data(:,Target_var_Pos));
Dum_var_Pos=cellfun(@(x) sum(strcmp(Dum_var_name,x))>0,Var_names,'UniformOutput',true);
Dummy=double(Data(:,Dum_var_Pos));
Exp_Var=double(Data(:,~(Target_var_Pos+Dum_var_Pos)));


%-------------------------Sample Period------------------------------------
Date=Data.Date(2:end);
Periodicity=4;

fyds=fix(Date(1));    %--- First Year of Data Set
if fyds>1300
    fyds=fyds-1300;
end
fqds=round(10*(Date(1)-fix(Date(1))));     %--- First quarter of Data Set
%lyds=92;   %--- Last Year of Data Set
%lqds=3;    %--- Last Quarter of Date Set
%-------------------------Forecasting Period for Evaluating Models---------
fyf=90;     %--- First Year of Forecasting
fqf=1;      %--- First Quarter of Forecast
lyf=95;     %--- Last Year of Forecast
lqf=1;      %--- Last Quarter of Forecast

%---------------------Preaparing Data--------------------------------------
T=[Target(1,:), Exp_Var(1,:)];
%T is transformation row vector for data
%so that [0 1 2 3 4 5]=[no change, Ln, Diff, Double Diff, Diff_Ln, Double Diff_Ln]
% Date=Target(2:end,1);

Target=Target(2:end,:);
Exp_Var=Exp_Var(2:end,:);
% F1=F1(2:end,2:end);
% F2=F2(2:end,2:end);
Dummy=Dummy(2:end,:);

[Target, Exp_Var, Dummy,Date]=Real_time(Target, Exp_Var, Dummy,Date);

[Target_Level_X12, Exp_Var]=deseasonal(Target, Exp_Var,Date,Periodicity);

[Target, Exp_Var, Dummy,Date]=transformation(Target_Level_X12, Exp_Var, Dummy, T,Date);

%--------------------------------------------------------------------------
dl=0;%length(Date)-length(Target);  %--- Real Time and transformation %% idon know what is this
sd=4*(fyf-fyds)+(fqf-fqds)-dl;     %--- Start Date to Forecasting
nd=4*(lyf-fyf)+(lqf-fqf);          %--- End Date to Forecasting
%--------------------------------------------------------------------------
model=nan(nd,horizon,300); % Each model Forcats+horizon
Desc_Model=cell(nd,300);%+horizon
%% find best Model
for i=1:nd%+horizon
    Model_Count=1;
    tic
    %     if OS==1 && i==nd+horizon
    %         OneStep=onestep;
    %     end
    %--------------------------Univariate Models---------------------------
    dum=Dummy(1:sd+i-1,:);
    %dum=zeros(size(Exp_Var,1),1);
    Target_=Target(1:sd-horizon-1+i,1);
    model(i,:,Model_Count)=ARdirect(Target_,dum,OneStep,horizon);   Desc_Model(i,Model_Count)={['ARdirect, Date, Upto, ' num2str(Date(sd+i))]};         Model_Count=Model_Count+1;
    model(i,:,Model_Count)=TAR(Target_,OneStep,horizon);           Desc_Model(i,Model_Count)={['TAR, Date, Upto, ' num2str(Date(sd+i))]};         Model_Count=Model_Count+1;
    model(i,:,Model_Count)=UM(Target_,OneStep,horizon);            Desc_Model(i,Model_Count)={['UM , Date, Upto, ' num2str(Date(sd+i))]};         Model_Count=Model_Count+1;
    model(i,:,Model_Count)=PureRW(Target_,OneStep,horizon);        Desc_Model(i,Model_Count)={['PureRW, Date, Upto, ' num2str(Date(sd+i))]};         Model_Count=Model_Count+1;
    model(i,:,Model_Count)=RWDrift(Target_,OneStep,horizon);        Desc_Model(i,Model_Count)={['RWDrift, Date, Upto, ' num2str(Date(sd+i))]};         Model_Count=Model_Count+1;
    model(i,:,Model_Count)=RWAO(Target_,OneStep,horizon);          Desc_Model(i,Model_Count)={['RWAO, Date, Upto, ' num2str(Date(sd+i))]};         Model_Count=Model_Count+1;
    %-----------------------Multivariate Models----------------------------
    %-------------------------------ARDL-----------------------------------
    Exp_Var_=Exp_Var(1:sd-horizon-1+i,:);
    
    %     [~,~, ~, ~, f_F1]=factoran(F1(1:sd-horizon-1+i,:), 3);
    %     [~,~, ~, ~, f_F2]=factoran(F2(1:sd-horizon-1+i,:), 2);
    
    %     Exp_Var_=[Exp_Var_];
    [~, CC]=size(Exp_Var_);
    %     Model_Count=6; % Counter of models
    for j=1:CC
        model(i,:,Model_Count)=ARDL2(Target_,Exp_Var_(:,j),dum,OneStep, horizon);
        Desc_Model(i,Model_Count)={['ARDL, Date, Upto, ' num2str(Date(sd+i)) ', Exo: ' num2str(j)]};
        Model_Count=Model_Count+1;
    end
    
    %-------------------------------VAR------------------------------------
    dataset=[Target_, Exp_Var_];
    
    [R,C]=size(dataset);
    nv=3;                 %%% number of variable in the model %%%
    v=nchoosek(2:C,nv);
    
    for k=1:size(v,1);
        
        Yraw=[dataset(:,1), dataset(:,v(k,:))];
        %-----------------------------------VAR----------------------------
        model(i,:,Model_Count)=BVAR(Yraw,dum,1,OneStep, horizon);     
        Desc_Model(i,Model_Count)={['BVAR, Date, Upto, ' num2str(Date(sd+i)) ', Exo: ' char(num2str(v(k,:)))]};         
        Model_Count=Model_Count+1;
        %----------------------------------TVP-VAR-------------------------
        
        model(i,:,Model_Count)=TVPVAR(Yraw,1,OneStep, horizon);
        Desc_Model(i,Model_Count)={['TVPVAR, Date, Upto, ' num2str(Date(sd+i)) ', Exo: ' char(num2str(v(k,:)))]};  
        Model_Count=Model_Count+1;
        %---------------------------------ARXs-----------------------------
        
        model(i,:,Model_Count)=ARXs(Yraw(:,1),Yraw(:,2:end),dum,OneStep, horizon);  
        Desc_Model(i,Model_Count)={['TVPVAR, Date, Upto, ' num2str(Date(sd+i)) ', Exo: ' char(num2str(v(k,:)))]}; 
        %Model_Count=Model_Count+1;
        %---------------------------------ARXd-----------------------------
        %
        %         model(i,:,Model_Count)=ARXd(dataset(:,1),dataset(:,v(k,:)),dum,4);
        %---------------------------------ARMAXs---------------------------
        %
        %         model(i,:,Model_Count)=ARMAXs(dataset(:,1),dataset(:,v(k,:)),dum,4,4);
    end
    %----------------------------------------------------------------------
    
    %     disp(i)
    tim=toc;
    %     round(tim*(nd+horizon-i));
    if round(rem(tim,1),3)==0
        disp(['Time Elapsed to compelete is: ' num2str(round(tim)) ' Second(s)']); %num2str(round(tim*(nd+horizon-i)))
    end
end

%Remove  Extra model
model(:,:,Model_Count+1:end)=[];

number_models=size(model,3);
%%
%**************************************************************************
%**************************************************************************
psedu_outofsample=nan(nd-horizon,horizon,number_models);
for i=1:horizon
    psedu_outofsample(:,i,:)=model(i:end-horizon+i-1,i,:);
end

outofsample=nan(horizon,number_models);
for i=1:horizon
    outofsample(i,:)=model(end,i,:);
end

%----------------Calculate RMSFE for all Models----------------------------
actual=Target(end-nd+horizon:end-1,1);
RMSFE1=nan(horizon,number_models);
for i=1:horizon
    for j=1:number_models
        RMSFE1(i,j)=RMSFE(actual, psedu_outofsample(:,i,j));
    end
end
% Sum of RMSE for each model over 4 horizon
RMSFE1_Total=sum(RMSFE1);
%[min_RMSFE, bestmodel0]=min(RMSFE1,[],2);
[min_RMSFE, bestmodel]=min(RMSFE1_Total);

disp(['the Best Model Base On RMSEF of psedu out of sample is ' Desc_Model{1,bestmodel}(1:min(strfind(Desc_Model{1,bestmodel},',')))]);
% for i=1:horizon
%     forecast_bestmodel(i,1)=outofsample(i,bestmodel(i,1));
% end


%---------------------Foracast Combination(Simple average)-----------------
for i=1:horizon
    l(i,1)=1;
    for k=1:number_models
        if RMSFE1(i,k)<=RMSFE1(i,1)
            psedu_outofsample_comb(:,i,l(i,1))=psedu_outofsample(:,i,k);
            outofsample_comb(i,l(i,1))=outofsample(i,k);
            l(i,1)=l(i,1)+1;
        end
    end
end

for i=1:horizon
    psedu_outofsample_combination(:,i)=mean(psedu_outofsample_comb(:,i,1:l(i,1)-1),3);
    outofsample_combination(i,1)=mean(outofsample_comb(i,1:l(i,1)-1),2);
end


for i=1:horizon
    RMSFE_comb(i,1)=RMSFE(actual, psedu_outofsample_combination(:,i));
end

% Output Decription
% actual: Observed Data
% model: (nd+horizon,horizon,300); % Each model Forcats
% Desc_Model:cell(nd+horizon,300); % Describe  eache model
% RMSFE1:nan(horizon,number_models);% RMSE of Eache model
% psedu_outofsample_combination: % best Out of sample Forcast
% outofsample_combination:
% Date: Date
Output=nan(nd+horizon,number_models*nd);
r=0;
for i=1:size(model,3)
    for j=1:size(model,1)
        r=r+1;
        Output(j:j+horizon-1,r)=squeeze(model(j,:,i));
    end
end

actual=[Target(sd+1:sd+nd);nan(horizon,1)];
Actual_date=[Date(sd+1:sd+nd);nan(horizon,1)];
for i=1:horizon
     Actual_date(nd+i)=Actual_date(nd+i-1)+0.01;
     if Actual_date(nd+i)-fix(Actual_date(nd+i))>0.115
         Actual_date(nd+i)=fix(Actual_date(nd+i))+1;
     end
end
% mode of Forcast
%M0 = mode(Output,2);
M1 = mean(Output,2,'omitnan');
M2=[nan(nd,1);forecast_bestmodel];
Output1=[Actual_date,actual,M1];
%----------------Writting forecast of best and combined model to plot fanchart---------
%{
xlswrite('Result.xlsx', 100*forecast_bestmodel, 'Bestmodel', 'B12:B15');
xlswrite('Result.xlsx', 100*min_RMSFE, 'Bestmodel', 'D12:D15');

xlswrite('Result.xlsx', 100*outofsample_combination, 'Combination', 'B7:B10');
xlswrite('Result.xlsx', 100*RMSFE_comb, 'Combination', 'D7:D10');

xlswrite('Result.xlsx', Date(end-13:end), 'Calculation', 'a2:a15');
xlswrite('Result.xlsx', 100*Target(end-9:end), 'Calculation', 'b2:b11');
xlswrite('Result.xlsx', Target_Level_X12(end-9:end), 'Calculation', 'c2:c11');

%----------------Writting forecast of best and combined model to plot fanchart---------

%}
save