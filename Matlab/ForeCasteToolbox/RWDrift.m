function [forecast]=RWDrift(inflation,one_step,horizon)
% it says inflation changes around the last observation with a drift
forecast=nan(1,horizon);
deltay=diff(inflation,1);
drift=(1/length(deltay))*sum(deltay);
for j=1:horizon
    forecast(1,j)=inflation(end)+j*drift;
end
if isnan(one_step)~=1
forecast(1,1)=one_step;
end
end

