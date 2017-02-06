function [forecast]=RWAO(inflation, one_step,horizon)
% it says inflation changes around the mean of last observations

y=mean(inflation(end-horizon+1:end,1));

% onestep=y;
% twostep=y;
% threestep=y;
% fourstep=y;
forecast=y*ones(1,horizon);%[onestep, twostep, threestep, fourstep];

if isnan(one_step)~=1
    forecast(1,1)=one_step;
end
end