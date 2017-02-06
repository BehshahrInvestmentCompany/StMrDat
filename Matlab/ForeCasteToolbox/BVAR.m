function [forecast]=BVAR(Yraw,dum, p,one_step,horizon)

%--------------------------------------------------------------------------
% Bayesian estimation, prediction and impulse response analysis in VAR
% models. Dependent on your choice of forecasting, the VAR model is:
%
% Iterated forecasts:
%     Y(t) = A0 + Y(t-1) x A1 + ... + Y(t-p) x Ap + e(t)
%
% so that in this case there are p lags of Y (from 1 to p).
%
% Direct h-step ahead foreacsts:
%     Y(t+h) = A0 + Y(t) x A1 + ... + Y(t-p+1) x Ap + e(t+h)
%
% so that in this case there are also p lags of Y (from 0 to p-1).
%
% In any of the two cases, the model is written as:
%
%                   Y(t) = X(t) x A + e(t)
%
% where e(t) ~ N(0,SIGMA), and A summarizes all parameters. Note that we
% also use the vector a which is defined as a=vec(A). 
%
%--------------------------------------------------------------------------


%------------------------------LOAD DATA-----------------------------------
% Load Quarterly US data on inflation, unemployment and interest rate, 
% 1953:Q1 - 2006:Q3
%load Yraw.dat;


% In any case, name the data you load 'Yraw', in order to avoid changing the
% rest of the code. Note that 'Yraw' is a matrix with T rows by M columns,
% where T is the number of time series observations (usually months or
% quarters), while M is the number of VAR dependent macro variables.
%----------------------------PRELIMINARIES---------------------------------
% Define specification of the VAR model
constant = 1;        % 1: if you desire intercepts, 0: otherwise 
%p = 2;               % Number of lags on dependent variables
forecasting = 1;     % 1: Compute h-step ahead predictions, 0: no prediction
forecast_method =1; % 0: Direct forecasts 
                     % 1: Iterated forecasts
h = horizon;               % Number of forecast periods

% Set prior for BVAR model:
prior = 1;  % prior = 1 --> Noninformative Prior
            % prior = 2 --> Minnesota Prior
            % prior = 3 --> Natural conjugate Prior
  
%--------------------------DATA HANDLING-----------------------------------
% Get initial dimensions of dependent variable
[Traw M] = size(Yraw);
%%%%Entering Dummy Variable%%%%
if sum(sum(dum(1:end-4,:)))==0
    d=0;
else
    d=1;
    s=sum(dum(1:end-4,:),1);
    m=1;
    for i=1:size(s,2)
        if s(1,i)~=0
            cn(1,m)=i;
            m=m+1;
        end
        
    end
end

% The model specification is different when implementing direct forecasts,
% compared to the specification when computing iterated forecasts.
if forecasting==1
    if h<=0
        error('You have set forecasting, but the forecast horizon choice is wrong')
    end    

    % Now create VAR specification according to forecast method
    if forecast_method==0       % Direct forecasts
        Y1 = Yraw(h+1:end,:);
        Y2 = Yraw(2:end-h,:);
        Traw = Traw - h - 1;
    elseif forecast_method==1   % Iterated forecasts
        Y1 = Yraw;
        Y2 = Yraw;
    else
        error('Wrong choice of forecast_method')
    end
else
   Y1 = Yraw;
   Y2 = Yraw;
end
        
% Generate lagged Y matrix. This will be part of the X matrix
Ylag = mlag2(Y2,p); % Y is [T x M]. ylag is [T x (Mp)]

% Now define matrix X which has all the R.H.S. variables (constant, lags of
% the dependent variable and exogenous regressors/dummies)
if constant
    switch d
        case 0
            X1 = [ones(Traw-p,1) Ylag(p+1:Traw,:)];
        case 1
            X1 = [ones(Traw-p,1) Ylag(p+1:Traw,:) dum(p+1:end-4,cn)];
    end
else
    switch d
        case 0
            X1 = Ylag(p+1:Traw,:);  %#ok<UNRCH>
        case 1
            X1 =[Ylag(p+1:Traw,:) dum(p+1:end-4,cn)]
    end
end

% Get size of final matrix X
[Traw3 K] = size(X1);

% Create the block diagonal matrix Z
Z1 = kron(eye(M),X1);

% Form Y matrix accordingly
% Delete first "LAGS" rows to match the dimensions of X matrix
Y1 = Y1(p+1:Traw,:); % This is the final Y matrix used for the VAR

% Traw was the dimesnion of the initial data. T is the number of actual 
% time series observations of Y and X
T = Traw - p;

%========= FORECASTING SET-UP:
% Now keep also the last "h" or 1 observations to evaluate (pseudo-)forecasts
%if forecasting==1
    %if forecast_method==0  % Direct forecasts, we only need to keep the 
        %Y = Y1(1:end-1,:);                             % last observation
        %X = X1(1:end-1,:);
        %Z = kron(eye(M),X);
        %T = T - 1;
    %else              % Iterated forecasts, we keep the last h observations
        %Y = Y1(1:end-h,:);
        %X = X1(1:end-h,:);
        %Z = kron(eye(M),X);
        %T = T - h;
    %end
%else
    Y = Y1;
    X = X1;
    Z = Z1;
%end

%--------------------------------PRIORS------------------------------------
% First get Ordinary Least Squares (OLS) estimators
A_OLS = inv(X'*X)*(X'*Y); % This is the matrix of regression coefficients
a_OLS = A_OLS(:);         % This is the vector of coefficients, i.e. it holds
                          % that a_OLS = vec(A_OLS)
SSE = (Y - X*A_OLS)'*(Y - X*A_OLS);
SIGMA_OLS = SSE./(T-K);

%-----------------Prior hyperparameters for bvar model
% Define hyperparameters
if prior == 1 % Noninformtive
    % I guess there is nothing to specify in this case!
    % Posteriors depend on OLS quantities
elseif prior == 2 % Minnesota
    A_prior = 0*ones(K,M);   
    a_prior = A_prior(:);
    
    % Hyperparameters on the Minnesota variance of alpha
    a_bar_1 = 0.7;
    a_bar_2 = 0.5;
    a_bar_3 = 10^2;
    
    % Now get residual variances of univariate p-lag autoregressions. Here
    % we just run the AR(p) model on each equation, ignoring the constant
    % and exogenous variables (if they have been specified for the original
    % VAR model)
    sigma_sq = zeros(M,1); % vector to store residual variances
    for i = 1:M
        % Create lags of dependent variable in i-th equation
        Ylag_i = mlag2(Yraw(:,i),p);
        Ylag_i = Ylag_i(p+1:Traw,:);
        % Dependent variable in i-th equation
        Y_i = Yraw(p+1:Traw,i);
        % OLS estimates of i-th equation
        alpha_i = inv(Ylag_i'*Ylag_i)*(Ylag_i'*Y_i);
        sigma_sq(i,1) = (1./(T-p+1))*(Y_i - Ylag_i*alpha_i)'*(Y_i - Ylag_i*alpha_i);
    end
    
    % Now define prior hyperparameters.
    % Create an array of dimensions K x M, which will contain the K diagonal
    % elements of the covariance matrix, in each of the M equations.
    V_i = zeros(K,M);
    
    % index in each equation which are the own lags
    ind = zeros(M,p);
    for i=1:M
        ind(i,:) = constant+i:M:K;
    end
    for i = 1:M  % for each i-th equation
        for j = 1:K   % for each j-th RHS variable
            if constant==1
                if j==1
                    V_i(j,i) = a_bar_3*sigma_sq(i,1); % variance on constant                
                elseif find(j==ind(i,:))>0
                    V_i(j,i) = a_bar_1./(p^2); % variance on own lags           
                else
                    for kj=1:M
                        if find(j==ind(kj,:))>0
                            ll = kj;                   
                        end
                    end
                    V_i(j,i) = (a_bar_2*sigma_sq(i,1))./((p^2)*sigma_sq(ll,1));           
                end
            else
                if find(j==ind(i,:))>0
                    V_i(j,i) = a_bar_1./(p^2); % variance on own lags
                else
                    for kj=1:M
                        if find(j==ind(kj,:))>0
                            ll = kj;
                        end                        
                    end
                    V_i(j,i) = (a_bar_2*sigma_sq(i,1))./((p^2)*sigma_sq(ll,1));            
                end
            end
        end
    end
    
    % Now V is a diagonal matrix with diagonal elements the V_i
    V_prior = diag(V_i(:));  % this is the prior variance of the vector a  
    
    % SIGMA is equal to the OLS quantity
    SIGMA = SIGMA_OLS;
    
elseif prior == 3 % Normal-Wishart (nat conj)
    % Hyperparameters on a ~ N(a_prior, SIGMA x V_prior)
    A_prior = 0*ones(K,M);   
    a_prior = A_prior(:);
    V_prior = 10*eye(K);
    % Hyperparameters on inv(SIGMA) ~ W(v_prior,inv(S_prior))
    v_prior = M;
    S_prior = eye(M);
    inv_S_prior = inv(S_prior);
end
    
%============================ POSTERIORS ==================================
%==========================================================================
    
%--------- Posterior hyperparameters of ALPHA and SIGMA with Diffuse Prior
if prior == 1
    % Posterior of alpha|Data ~ Multi-T(kron(SSE,inv(X'X)),alpha_OLS,T-K)
    V_post = inv(X'*X);
    a_post = a_OLS;
    A_post = reshape(a_post,K,M);
    
    % posterior of SIGMA|Data ~ inv-Wishart(SSE,T-K)
    S_post = SSE;
    v_post = T-K;
    
    % Now get the mean and variance of the Multi-t marginal posterior of alpha
    alpha_mean = a_post;
    alpha_var = (1/(v_post - M - 1))*kron(S_post,V_post);
    
%--------- Posterior hyperparameters of ALPHA and SIGMA with Minnesota Prior
elseif prior == 2
    % ******Get all the required quantities for the posteriors       
    V_post = inv( inv(V_prior) + kron(inv(SIGMA),X'*X) );
    a_post = V_post * ( inv(V_prior)*a_prior + kron(inv(SIGMA),X'*X)*a_OLS );
    A_post = reshape(a_post,K,M);
     
    % In this case, the mean is a_post and the variance is V_post
   alpha_mean=a_post;
%--------- Posterior hyperparameters of ALPHA and SIGMA with Normal-Wishart Prior
elseif prior == 3
    % ******Get all the required quantities for the posteriors       
    % For alpha
    V_post = inv( inv(V_prior) + X'*X );
    A_post = V_post * ( inv(V_prior)*A_prior + X'*X*A_OLS );
    a_post = A_post(:);
    
    % For SIGMA
    S_post = SSE + S_prior + A_OLS'*X'*X*A_OLS + A_prior'*inv(V_prior)*A_prior - A_post'*( inv(V_prior) + X'*X )*A_post;
    v_post = T + v_prior;
    
    % Now get the mean and variance of the Multi-t marginal posterior of alpha
    alpha_mean = a_post;
    alpha_var = (1/(v_post - M - 1))*kron(S_post,V_post); 
end



%%% AIC %%%

%AIC=log(det(S_post))+(2*(M*(p*M+1)))/T;
%SIC=log(det(S_post))+(log(T))*(M*(p*M+1))/T;



%======================= Forecasting 1 step to 4steps =============================
%==========================================================================
%X_tplus1 = [1 Y(T,:) X(T,2:M*(p-1)+1)];

Xf=Yraw;
for k=1:4
    X_tplus = mlag3(Xf,p);
    switch d
        case 0
            X_tplus1=[1, X_tplus(end,:)];
        case 1
            X_tplus1=[1, X_tplus(end,:) dum(end+k-4,cn)];
    end
    forecast1(k,:) = X_tplus1*A_post;
    %-------------------------------------------%
    if k==1 && isnan(one_step)~=1
        forecast1(k,1)=one_step;
    end
    %-------------------------------------------%
    Xf=[Xf; forecast1(k,:)];

end

    forecast=forecast1(1:4,1)';

%-----------------------Implicit Functions-------------------------------

function [Xlag] = mlag2(X,p)
%MLAG2 Summary of this function goes here
%   Detailed explanation goes here
[Traw,N]=size(X);
Xlag=zeros(Traw,N*p);
for ii=1:p
    Xlag(p+1:Traw,(N*(ii-1)+1):N*ii)=X(p+1-ii:Traw-ii,1:N);
end


% %OR:
% [Traw,N]=size(X);
% Xlag=zeros(Traw,N,p);
% for ii=1:p
%     Xlag(p+1:Traw,:,ii)=X(p+1-ii:Traw-ii,:);
% end
% Xlag=Xlag(:,:);

function [Xlagf] = mlag3(X,p)
%MLAG2 Summary of this function goes here
%   Detailed explanation goes here
[Traw,N]=size(X);
Xlag=zeros(Traw,N*p);
for ii=1:p
    Xlag(p+1:Traw,(N*(ii-1)+1):N*ii)=X(p+1-ii:Traw-ii,1:N);
    Xlagf=[X(:,:), Xlag(:,1:N*(p-1))];
end


function A = wish(h,n)
% Command:  s=wish(h,n)
% Purpose:  Draws an m x m matrix from a wishart distribution
%           with scale matrix h and degrees of freedom nu = n.
%           This procedure uses Bartlett's decomposition.
% Inputs:   h     -- m x m scale matrix.
%           n     -- scalar degrees of freedom.
% Outputs:  s     -- m x m matrix draw from the wishart
%                    distribution.
% Note: Parameterized so that mean is n*h

A = chol(h)'*randn(size(h,1),n);
A = A*A';

function [y] = bvardgp(T,N,L,PHI,PSI)
%--------------------------------------------------------------------------
%   PURPOSE:
%      Get matrix of Y generated from a VAR model
%--------------------------------------------------------------------------
%   INPUTS:
%     T     - Number of observations (rows of Y)
%     N     - Number of series (columns of Y)
%     L     - Number of lags
%
%   OUTPUT:
%     y     - [T x N] matrix generated from VAR(L) model
% -------------------------------------------------------------------------

randn('seed',sum(100*clock));
rand('seed',sum(100*clock));
%-----------------------PRELIMINARIES--------------------
if nargin==0;
    T = 400;           %Number of time series observations (T)
    N = 6;             %Number of cross-sectional observations (N)
    L = 1;             %Lag order

    PHI = [1 1 1 1 1 1;...
          .5 0 0 0 0 0;...
           0 .5 0 0 0 0;...
           0 0 .5 0 0 0;...
           0 0 0 .5 0 0;...
           0 0 0 0 .5 0;...
           0 0 0 0 0 .5];
    
PSI = [1.70    0.02    0.16    0.08    0.19   -0.03;...
       0.02    0.66    0.16    0.03    0.26    0.07;...
       0.16    0.16    0.94    0.25    0.06    0.24;...
       0.08    0.03    0.25    1.45    0.06    0.38;...
       0.19    0.26    0.06    0.06    0.64    0.16;...
      -0.03    0.07    0.24    0.38    0.16    1.33;];
end

%---------------------------------------
% Ask user if a constant is desired
f = input('Do you want to generate a VAR model with a constant? <y/n>: ','s');
if strcmp(f,'yes') || strcmp(f,'y') % compare strings. If f = 'yes' then
    const = 1;
    disp(['VAR with ' num2str(L) '-lag(s) and an intercept generated'])
else                % elseif f ='no' then
    const=0;
    disp(['VAR with ' num2str(L) '-lag(s) with NO intercept generated'])
end

%----------------------GENERATE--------------------------
% Set storage in memory for y
% First L rows are created randomly and are used as 
% starting (initial) values 
y =[rand(L,N) ; zeros(T-L,N)];

% Now generate Y from VAR (L,PHI,PSI)
for nn = L:T
    u = chol(PSI)'*randn(N,1);
    ylag = mlag2(y,L);
    y(nn,:) = [const ylag(nn,:)]*PHI + u';
end
% %-------------------------MLE----------------------------
% % Now we can estimate PHI,PSI using OLS. First transform
% % data to get Y and X matrices. X = constant + lagged(Y)
% ylag = mlag(y,L); % Create lagged Y matrix 
% 
% if const==1  % Create X matrix with or without a constant
%     xmat = [ones(T-L,1) ylag(L+1:T,:)];
% else
%     xmat = ylag(L+1:T,:);
% end
% % Chop off first L obs from Y, to match dimensions
% % of X matrix
% ymat = y(L+1:T,:); 
% 
% % Now get MLE quantities.
% PHI_M = inv(xmat'*xmat)*(xmat'*ymat);
% SSE = (ymat - xmat*PHI_M)'*(ymat - xmat*PHI_M);

