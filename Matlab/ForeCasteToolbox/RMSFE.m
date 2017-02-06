function [RMSFE,ForecastError]=RMSFE(actual, forecast)
% Root Mean Square Forecaste Error
T=length(actual);
actual=actual(:);
forecast=forecast(:);
e=actual-forecast;
e2=(e'*e)/T;
ForecastError=e(:);
RMSFE=(e2)^0.5;