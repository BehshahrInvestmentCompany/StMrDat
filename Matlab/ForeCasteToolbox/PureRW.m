function [forecast]=PureRW(inflation,one_step,horizon)
% it says inflation changes around the last observation 
% onestep=inflation(end,1);
% twostep=inflation(end,1);
% threestep=inflation(end,1);
% fourstep=inflation(end,1);

forecast=inflation(end,1)*ones(1,horizon);%[onestep, twostep, threestep, fourstep];
if isnan(one_step)~=1
    forecast(1,1)=one_step;
end
end