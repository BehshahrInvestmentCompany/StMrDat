% [$-fa-IR,300]#,##0
%{
f=fullfile(cd,'\IRIS_Tbx');
addpath(f)
irisstartup
%}
% Version 2
% Base on Karami and Bayat
% Modified By P.Davoudi: Pedram.Davoudi@gmail.com
clear
clc
Irann=1;
JustAr=1;
%****************************
onestep=NaN; % Idont know it
OneStep=nan; % Idont know it
horizon=4; % Forecast Horizon
if Irann==1
    isPersianDate=1;
    % %---------------------Read Data from Excell file---------------------------
    Data=dataset('xls','Input\Data.xlsx','sheet','Data');
    Target_var_names={'CPI'};%:,'Gold','USD','SI','GDP'};
    Dum_var_name={'D1','D2'};
    nd=20; % Periods to find the best models
    %the first row is the transformations
    FRT=1;
    %%
else % Commodities
    isPersianDate=0;
    % %---------------------Read Data from Excell file---------------------------
    Data=dataset('xls','Input\Mar.xlsx');%,'sheet','Data');
    Target_var_names={'CU'};% 'oil' 'Al' 'Co' 'Or'};
    Dum_var_name={};%{'D1','D2'}; periodical=Data.Date(1); % 4 for quarterly data and 12 for monthly and 1 for anual
    nd=10; % Periods to find the best models
    %the first row is the transformations
    FRT=1;
end
for Tvarnameindx=1:length(Target_var_names)
    Target_var_name=Target_var_names{Tvarnameindx};
    % ---------------------------------------------------
    Var_names= Data.Properties.VarNames;
    Target_var_Pos=strcmp(Var_names,Target_var_name);
    Target=double(Data(:,Target_var_Pos));
    Date_var_Pos=strcmp(Var_names,'Date');
    
    Dum_var_Pos=cellfun(@(x) sum(strcmp(Dum_var_name,x))>0,Var_names,'UniformOutput',true);
    Dummy=double(Data(:,Dum_var_Pos));
    if all(isempty(Dummy))
        Dummy=zeros(size(Dummy,1),1);
    end
    Exp_Var=double(Data(:,~(Date_var_Pos+Target_var_Pos+Dum_var_Pos)));
    Exp_Var_names=Var_names(~(Date_var_Pos+Target_var_Pos+Dum_var_Pos));
    Date=double(Data(:,Date_var_Pos));
    
    % dataset=[ Exp_Var_];
    C=length(Exp_Var_names);
    nv=3;                 %%% number of variable in the model %%%
    v=nchoosek(1:C,nv);
    
    %---------------------Preaparing Data--------------------------------------
    T=[Target(1,:), Exp_Var(1,:)];
    if FRT==1
        Target=Target(2:end,:);
        Exp_Var=Exp_Var(2:end,:);
        % F1=F1(2:end,2:end);
        % F2=F2(2:end,2:end);
        periodicity=Date(1); % 4 for quarterly data and 12 for monthly and 1 for anual
        Dummy_Orig=Dummy(2:end,:);
        Date=Date(2:end);
    elseif FRT==0
        Dummy_Orig=Dummy(1:end,:);
        periodicity=12;
        T=410*ones(size(T));
    end
    if isPersianDate
        Date_=Cal_conv(Date,0);
    else
        Date_=Date;
    end
    %T is transformation row vector for data
    %so that [0 1 2 3 4 5]=[no change, Ln, Diff, Double Diff, Diff_Ln, Double Diff_Ln]
    Trans_mod={'no\_change', 'Ln', 'Diff', 'Double\_Diff', 'Diff\_Ln', 'Double\_Diff\_Ln'};
    Tt2=floor(log(T(1))./log(10))+1;
    if ~isfinite(Tt2)
        Tt2=1;
    end
    temp_=[];
    for tt=1:Tt2
        temp_=[temp_;(mod(T,10^(tt))-mod(T,10^(tt-1)))./10^(tt-1)];
    end
    T=temp_;
    clear temp_
    
    [Target, Exp_Var, Dummy_,Date_]=Real_time(Target, Exp_Var, Dummy_Orig,Date_);
    
    [Target_Level_X12, Exp_Var_B]=deseasonal(Target, Exp_Var,Date_,periodicity);
    Target_Historical=Target;
    
    %% Simulate for Model Selection
    model=nan(nd,horizon,Tt2,500); % Each model Forcats+horizon
    Desc_Model=cell(nd,Tt2,500);%+horizon
    Full_model=nan(horizon,tt,500); % Fulll Samples
    Fin_Desc_Model=cell(tt,500);
    
    for tt=1:Tt2
        %tt=2
        [Target, Exp_Var, Dummy,Date]=transformation(Target_Level_X12, Exp_Var_B, Dummy_, T(tt,:),Date_);
        %-------------------------Sample Period------------------------------------
        %--- First Year of Data Set
        fyds=str2double(datestr(Date(1),'yyyy'));
        %--- First Month of Data Set
        fqds=str2double(datestr(Date(1),'mm'));
        %--- Last Year of Historical Data
        lyf=str2double(datestr(Date(end),'yyyy'));
        %--- Last Month of Historical Data
        lqf=str2double(datestr(Date(end),'mm'));
        %--- First Year of Sampling
        fyf=str2double(datestr(Date(end)-nd*365/periodicity,'yyyy'));
        %--- First Month of Sampling
        fqf=str2double(datestr(Date(end)-nd*365/periodicity,'mm'));
        dl=0;%length(Date)-length(Target);  %--- Real Time and transformation %% idon know what is this
        sd=periodicity*((fyf-fyds)+(fqf-fqds)/12)-dl;     %--- Start Date to Forecasting
        %nd=4*(lyf-fyf)+(lqf-fqf);          %--- End Date to Forecasting
        %--------------------------------------------------------------------------
        %-------------------------Forecasting Period for Evaluating Models---------
        
        
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
            model(i,:,tt,Model_Count)=ARdirect(Target_,dum,OneStep,horizon);   Desc_Model(i,tt,Model_Count)={['ARdirect, Date, Upto, ' num2str(Date(sd+i)) ', Transfomation mod,' Trans_mod{tt}]};         Model_Count=Model_Count+1;
            model(i,:,tt,Model_Count)=TAR(Target_,OneStep,horizon);            Desc_Model(i,tt,Model_Count)={['TAR, Date, Upto, ' num2str(Date(sd+i)) ', Transfomation mod, ' Trans_mod{tt}]};              Model_Count=Model_Count+1;
            model(i,:,tt,Model_Count)=UM(Target_,OneStep,horizon);             Desc_Model(i,tt,Model_Count)={['UM , Date, Upto, ' num2str(Date(sd+i)) ', Transfomation mod, ' Trans_mod{tt}]};              Model_Count=Model_Count+1;
            model(i,:,tt,Model_Count)=PureRW(Target_,OneStep,horizon);         Desc_Model(i,tt,Model_Count)={['PureRW, Date, Upto, ' num2str(Date(sd+i)) ', Transfomation mod, ' Trans_mod{tt}]};           Model_Count=Model_Count+1;
            model(i,:,tt,Model_Count)=RWDrift(Target_,OneStep,horizon);        Desc_Model(i,tt,Model_Count)={['RWDrift, Date, Upto, ' num2str(Date(sd+i)) ', Transfomation mod, ' Trans_mod{tt}]};          Model_Count=Model_Count+1;
            model(i,:,tt,Model_Count)=RWAO(Target_,OneStep,horizon);           Desc_Model(i,tt,Model_Count)={['RWAO, Date, Upto, ' num2str(Date(sd+i)) ', Transfomation mod, ' Trans_mod{tt}]};             Model_Count=Model_Count+1;
            %-----------------------Multivariate Models----------------------------
            if JustAr==0
                %-------------------------------ARDL-----------------------------------
                Exp_Var_=Exp_Var(1:sd-horizon-1+i,:);
                
                %     [~,~, ~, ~, f_F1]=factoran(F1(1:sd-horizon-1+i,:), 3);
                %     [~,~, ~, ~, f_F2]=factoran(F2(1:sd-horizon-1+i,:), 2);
                
                %     Exp_Var_=[Exp_Var_];
                [~, CC]=size(Exp_Var_);
                %     Model_Count=6; % Counter of models
                for j=1:CC
                    model(i,:,tt,Model_Count)=ARDL2(Target_,Exp_Var_(:,j),dum,OneStep, horizon); Desc_Model(i,tt,Model_Count)={['ARDL, Date, Upto, ' num2str(Date(sd+i)) ', Transfomation mod, ' Trans_mod{tt} ', Exo: ' Exp_Var_names{j}]};         Model_Count=Model_Count+1;
                end
                
                %-------------------------------VAR------------------------------------
                %dataset=[Target_, Exp_Var_];
                
                
                for k=1:size(v,1);
                    
                    Yraw=[Target_, Exp_Var_(:,v(k,:))];
                    %-----------------------------------VAR----------------------------
                    model(i,:,tt,Model_Count)=BVAR(Yraw,dum,1,OneStep, horizon);
                    %        [ i,tt,Model_Count]
                    Desc_Model(i,tt,Model_Count)={['BVAR, Date, Upto, ' num2str(Date(sd+i)) ', Transfomation mod, ' Trans_mod{tt} ', Exo: ' Exp_Var_names{v(k,:)}]};
                    Model_Count=Model_Count+1;
                    %----------------------------------TVP-VAR-------------------------
                    
                    model(i,:,tt,Model_Count)=TVPVAR(Yraw,1,OneStep, horizon);                Desc_Model(i,tt,Model_Count)={['TVPVAR, Date, Upto, ' num2str(Date(sd+i)) ', Transfomation mod, ' Trans_mod{tt} ', Exo: ' Exp_Var_names{v(k,:)}]};         Model_Count=Model_Count+1;
                    %---------------------------------ARXs-----------------------------
                    
                    model(i,:,tt,Model_Count)=ARXs(Yraw(:,1),Yraw(:,2:end),dum,OneStep, horizon);     Desc_Model(i,tt,Model_Count)={['TVPVAR, Date, Upto, ' num2str(Date(sd+i)) ', Transfomation mod, ' Trans_mod{tt} ', Exo: ' Exp_Var_names{v(k,:)}]};         Model_Count=Model_Count+1;
                    %---------------------------------ARXd-----------------------------
                    %
                    %         model(i,:,tt,Model_Count)=ARXd(dataset(:,1),dataset(:,v(k,:)),dum,4);
                    %---------------------------------ARMAXs---------------------------
                    %
                    %         model(i,:,tt,Model_Count)=ARMAXs(dataset(:,1),dataset(:,v(k,:)),dum,4,4);
                end
                %----------------------------------------------------------------------
            end
            %     disp(i)
            tim=toc;
            %     round(tim*(nd+horizon-i));
            if round(rem(tim,1),3)==0
                disp(['Time Elapsed to compelete is: ' num2str(round(tim)) ' Second(s)']); %num2str(round(tim*(nd+horizon-i)))
            end
        end
        
        
        %% Simulate for full Sample data
        
        Model_Count=1;
        %     if OS==1 && i==nd+horizon
        %         OneStep=onestep;
        %     end
        %--------------------------Univariate Models---------------------------
        
        %dum=zeros(size(Exp_Var,1),1);
        %         Target_=Target(1:sd-horizon-1+i,1);
        Target_=Target;
        dum=Dummy_Orig(1:length(Target_)+horizon,:);
        Full_model(:,tt,Model_Count)=ARdirect(Target_,dum,OneStep,horizon);   Fin_Desc_Model(tt,Model_Count)={['ARdirect, Trans:' Trans_mod{tt}]}; Model_Count=Model_Count+1;
        Full_model(:,tt,Model_Count)=TAR(Target_,OneStep,horizon);            Fin_Desc_Model(tt,Model_Count)={['TAR, Trans:' Trans_mod{tt}]}; Model_Count=Model_Count+1;
        Full_model(:,tt,Model_Count)=UM(Target_,OneStep,horizon);            Fin_Desc_Model(tt,Model_Count)={['UM , Trans:' Trans_mod{tt}]};   Model_Count=Model_Count+1;
        Full_model(:,tt,Model_Count)=PureRW(Target_,OneStep,horizon);        Fin_Desc_Model(tt,Model_Count)={['PureRW, Trans:' Trans_mod{tt}]};   Model_Count=Model_Count+1;
        Full_model(:,tt,Model_Count)=RWDrift(Target_,OneStep,horizon);       Fin_Desc_Model(tt,Model_Count)={['RWDrift, Trans:' Trans_mod{tt}]};   Model_Count=Model_Count+1;
        Full_model(:,tt,Model_Count)=RWAO(Target_,OneStep,horizon);          Fin_Desc_Model(tt,Model_Count)={['RWAO, Trans:' Trans_mod{tt}]};   Model_Count=Model_Count+1;
        %-----------------------Multivariate Models----------------------------
        if JustAr==0
            %-------------------------------ARDL-----------------------------------
            %  Exp_Var_=Exp_Var(1:sd-horizon-1+i,:);
            Exp_Var_=Exp_Var;
            %     [~,~, ~, ~, f_F1]=factoran(F1(1:sd-horizon-1+i,:), 3);
            %     [~,~, ~, ~, f_F2]=factoran(F2(1:sd-horizon-1+i,:), 2);
            
            %     Exp_Var_=[Exp_Var_];
            [~, CC]=size(Exp_Var_);
            %     Model_Count=6; % Counter of models
            for j=1:CC
                Full_model(:,tt,Model_Count)=ARDL2(Target_,Exp_Var_(:,j),dum,OneStep, horizon); Fin_Desc_Model(tt,Model_Count)={['ARDL, Trans:' Trans_mod{tt} ', Exo: ' Exp_Var_names{j}]};   Model_Count=Model_Count+1;
            end
            
            %-------------------------------VAR------------------------------------
            %dataset=[ Exp_Var_]; %Target_,
            
            % [R,C]=size(dataset);
            %nv=3;                 %%% number of variable in the model %%%
            %v=nchoosek(1:C,nv);
            
            for k=1:size(v,1);
                
                Yraw=[Target_(:,1), Exp_Var_(:,v(k,:))];
                %-----------------------------------VAR----------------------------
                Full_model(:,tt,Model_Count)=BVAR(Yraw,dum,1,OneStep, horizon);       Fin_Desc_Model(tt,Model_Count)={['BVAR, Trans:' Trans_mod{tt} ', Exo: ' Exp_Var_names{v(k,:)}]};   Model_Count=Model_Count+1;
                %----------------------------------TVP-VAR-------------------------
                
                Full_model(:,tt,Model_Count)=TVPVAR(Yraw,1,OneStep, horizon);         Fin_Desc_Model(tt,Model_Count)={['TVPVAR, Trans:' Trans_mod{tt} ', Exo: ' Exp_Var_names{v(k,:)}]};   Model_Count=Model_Count+1;
                %---------------------------------ARXs-----------------------------
                
                Full_model(:,tt,Model_Count)=ARXs(Yraw(:,1),Yraw(:,2:end),dum,OneStep, horizon);       Fin_Desc_Model(tt,Model_Count)={['TVPVAR, Trans:' Trans_mod{tt} ', Exo: ' Exp_Var_names{v(k,:)}]};   Model_Count=Model_Count+1;
                %---------------------------------ARXd-----------------------------
                %
                %         model(i,:,tt,Model_Count)=ARXd(dataset(:,1),dataset(:,v(k,:)),dum,4);
                %---------------------------------ARMAXs---------------------------
                %
                %         model(i,:,tt,Model_Count)=ARMAXs(dataset(:,1),dataset(:,v(k,:)),dum,4,4);
            end
            %----------------------------------------------------------------------
        end
    end
    %%
    %  Fin_Desc_Model(:,1,:)=[]; % remove extra
    model_Backup=model;
    Full_model_Backup=Full_model;
    model(:,:,:,Model_Count:end)=[];%squeeze(isnan(model(1,1,:)))
    Full_model(:,:,Model_Count:end)=[];
    Desc_Model(:,:,Model_Count:end)=[];
    Fin_Desc_Model(:,Model_Count:end)=[];
    %%
    %Remove  Extra model
    
    Model_Count=Model_Count-1;
    
    % Detrandorm forecasts
    for i=1:nd
        for j=1:Model_Count
            for tt=1:Tt2
                base=Target_Level_X12(sd-horizon-1+i,1);
                model(i,:,tt,j)=DeTransformation(squeeze(model(i,:,tt,j)),base, T(tt,1));
            end
        end
    end
    
    base=Target_Level_X12(end,1);
    for j=1:Model_Count
        for tt=1:Tt2
            Full_model(:,tt,j)=DeTransformation(squeeze(Full_model(:,tt,j)),base, T(tt,1));
        end
    end
    %%
    %**************************************************************************
    %**************************************************************************
    % psedu_outofsample=nan(nd-horizon,horizon,Tt2,number_models);
    psedu_outofsample=nan(nd,horizon,Tt2,Model_Count);
    
    
    for i=1:horizon
        psedu_outofsample(:,i,:,:)=model(:,i,:,:);%model(i:end-horizon+i-1,i,:,:);
    end
    
    %     outofsample=nan(horizon,Tt2,Model_Count);
    %     for i=1:horizon
    %         outofsample(i,:,:)=model(end,i,:,:);
    %     end
    
    %----------------Calculate RMSFE for all Models----------------------------
    %%
    RMSFE1=nan(nd,Tt2,Model_Count);
    forecast_Errors=[];
    for tt=1:Tt2
        for i=1:nd
            for j=1:Model_Count
                %[tt i j]
                % Target_=Target(1:sd-horizon-1+i,1);
                actual=Target_Level_X12(sd+i-horizon-1+1:sd-horizon-1+i+4,1);
                [RMSFE1(i,tt,j), forecast_Error]=RMSFE(actual, squeeze(psedu_outofsample(i,:,tt,j)));
                forecast_Errors=[forecast_Errors;forecast_Error,tt*ones(size(forecast_Error))];
            end
        end
    end
    % Sum of RMSE for each model over 4 horizon
    RMSFE1_Total=squeeze(sum(RMSFE1,1));
    min_RMSFE=min(min(RMSFE1_Total));
    [R,C]=find(RMSFE1_Total==min_RMSFE);
    [~, bestmodel0]=min(RMSFE1,[],2);
    
    disp(['the Best Model to Forecaste ' Target_var_name ' Base On RMSEF of psedu out of sample is ' Desc_Model{1,R,C}] );% Desc_Model{1,bestmodel}(1:min(strfind(Desc_Model{1,bestmodel},',')))]);
    %     for i=1:horizon
    %         forecast_bestmodel(i,1)=outofsample(i,bestmodel0(i,1));
    %     end
    %%
    % Output Decription
    % actual: Observed Data
    % model: (nd+horizon,horizon,300); % Each model Forcats for model
    % Comparison
    % Full_model: Each model Forcats for Ful Sample
    % Desc_Model:cell(nd+horizon,300); % Describe  eache model
    % RMSFE1:nan(horizon,Model_Count);% RMSE of Eache model
    % Fin_Desc_Model:cell(nd+horizon,300); % Describe  eache model
    % psedu_outofsample_combination: % best Out of sample Forcast
    % outofsample_combination:
    % Date: Date
    Historical_Number=5;
    R_Full_model=reshape(Full_model,horizon,[]);
    %Result_Database=mat2dataset([repmat(Target_Level_X12(end-Historical_Number+1:end).',size(R_Full_model,2),1),R_Full_model.',RMSFE1_Total(:)],'varnames',{'L1', 'L2', 'L3' , 'L4', 'L5','H1', 'H2', 'H3' , 'H4', 'RMSEF' });
    Result_Database=mat2dataset([repmat(Target_Historical(end-Historical_Number+1:end).',size(R_Full_model,2),1),R_Full_model.',RMSFE1_Total(:)],...
        'varnames',[ mat2cell(datestr(Date_(end-Historical_Number+1:end)),ones(5,1),11).',{'H1', 'H2', 'H3' , 'H4', 'RMSEF' }]);
    
    Result_Database.Descritption=Fin_Desc_Model(:);
    Result_Database = sortrows(Result_Database,'RMSEF','ascend');
    
    Top_count=5;
    figure;
    plot(R_Full_model,'Color',[0.9,0.9,0.9])
    if fix(Model_Count/10)>5
        Fin_fore=double(Result_Database(1:fix(Model_Count/10),1:Historical_Number+horizon));
    else
        Fin_fore=double(Result_Database(1:fix(Model_Count),1:Historical_Number+horizon));
    end
    plot(Fin_fore(1:Top_count,:).')
    hold on
    plot(mode(Fin_fore))
    plot(mean(Fin_fore))
    leg=cell(1,Top_count+2);
    %     for ll=1:Top_count
    %         leg{1,ll}=['#' num2str(ll)];
    %     end
    leg(1:Top_count)=Result_Database.Descritption(1:Top_count);
    leg{1,Top_count+1}='Mode';
    leg{1,Top_count+2}='Mean';
    legend (leg,'Location','northwest');
    title(['4 Horizon Forecaste of ' Target_var_name ' __ From ' datestr(Date(end))]);
    JmP=Date(end)-Date(end-1);
    
    ax=gca;
    ax.XTickLabel =(datestr(Date(end-4):JmP:Date(end-4)+8*JmP,'yyyy-mmm'));
    hold off
    if ~exist('Output','dir')
        mkdir('Output')
    end
    saveas(gcf,['Output\Lev_' Target_var_name],'bmp');
    % close gcf
    % Growth of target variable
    figure;
    %  plot(diff(log(R_Full_model)),'Color',[0.9,0.9,0.9])
    D_Fin_fore=diff(log(Fin_fore),1,2);%double(Result_Database(1:fix(Model_Count/10),1:horizon));
    plot(D_Fin_fore(1:Top_count,:).')
    hold on
    plot(mode(D_Fin_fore))
    plot(mean(D_Fin_fore))
    leg=cell(1,Top_count+2);
    leg(1:Top_count)=Result_Database.Descritption(1:Top_count);
    leg{1,Top_count+1}='Mode';
    leg{1,Top_count+2}='Mean';
    legend (leg,'Location','northwest');
    title(['4 Horizon Forecaste of ' Target_var_name ' Growth __ From ' datestr(Date(end),'yy-mmm')]);
    JmP=Date(end)-Date(end-1);
    
    ax=gca;
    ax.XTickLabel =(datestr(Date(end-4):JmP:Date(end-4)+8*JmP,'yyyy-mmm'));
    hold off
    saveas(gcf,['Output\Grt_' Target_var_name],'bmp');
    % close gcf
    
    % forecast_Errors=sort(forecast_Errors);
    % ksdensity(forecast_Errors(forecast_Errors(500:end-500,2)==1,1))
    % ksdensity(forecast_Errors(500:end-500))
    % title('Forecaste Error Distribution');
    save(['Output\' Target_var_name '.mat']);
    export(Result_Database,'xlsfile',['Output\' Target_var_name '.xlsx']);
end