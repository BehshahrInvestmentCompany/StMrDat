function [forecast]=UM(inflation, one_step,horizon)
% it says inflation changes around a constatnt term with white noise
ucmean=mean(inflation);

forecast=ucmean*ones(1,horizon);%[ucmean,ucmean,ucmean,ucmean];
if isnan(one_step)~=1
    forecast(1,1)=one_step;
end
end