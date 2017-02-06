function [forecast]=TVPVAR(Yraw,k, one_step,horizon)



%we have state space representation as:

%y(t)=x(t)*beta(t)+w(t)   w~N(0,Q);
%beta(t+1)=F*beta(t)+v(t+1)             v~N(0,R) ;


%T:number of observation; N=number of variabales (N=4); y:T*N
% X:T*N;  beta:(N+1)*(N+1); w:T*N ; Q:N*N; v:(N+1)*(N+1);
% R:N*N; F:(N+1)*(N+1)

%K is the number of lags in each equation of VAR.


%-------------------------------- Constructing VAR -----------------------------------
data3=Yraw;

[T,N]=size(data3);

y=[];
x=[];
for m=1:N
    y(:,m)=data3(k+1:T,m);
    for j=1:k
        x(:,k*(m-1)+j)=data3(k+1-j:T-j,m);
    end
end
[u,v]=size(x);
x=[ones(u,1) x];

%--------------------------------Parameter Values from OLS in the First Equation of VAR--------------------------------------

beta(:,:,1)=zeros(5,4);
betat(:,:,1)=zeros(5,4);

yy=data3(:,1);
xx=lagmatrix([ones(length(yy),1) data3],1);
yy=yy(2:end,:);
xx=xx(2:end,:);
b=inv(xx'*xx)*xx'*yy;
resid1=yy-xx*b;
sigma1=((length(yy)-length(b))/length(yy))*(resid1'*resid1);
varb1=sigma1*inv(xx'*xx);

beta(1,1,1)=b(1,1); beta(2,1,1)=b(2,1); beta(3,1,1)=b(3,1);  beta(4,1,1)=b(4,1);  beta(5,1,1)=b(5,1);

p1(:,:,1)=varb1;
Q=sigma1;
R=varb1;

F1=eye(5);
%--------------------------------kalman Filter for Parameter Estimation in the First Equation of VAR--------------------------------------
for t=1:T-1
    p1(:,:,t+1)=F1*(p1(:,:,t)-(p1(:,:,t)*x(t,:)'*inv(x(t,:)*p1(:,:,t)*x(t,:)'+Q)*x(t,:)*p1(:,:,t)))*F1'+R*R';
    beta(:,1,t+1)=F1*beta(:,1,t)+F1*p1(:,:,t)*x(t,:)'*inv(x(t,:)*p1(:,:,t)*x(t,:)'+Q)*(y(t,1)-x(t,:)*beta(:,1,t));
    betat(:,1,t)=beta(:,1,t)+p1(:,:,t)*x(t,:)'*inv(x(t,:)*p1(:,:,t)*x(t,:)'+Q)*(y(t,1)-x(t,:)*beta(:,1,t));
end

%--------------------------------Parameter Values from OLS in the second Equation of VAR--------------------------------------
yy=data3(:,2);
xx=lagmatrix([ones(length(yy),1) data3],1);
yy=yy(2:end,:);
xx=xx(2:end,:);
b=inv(xx'*xx)*xx'*yy;
resid2=yy-xx*b;
sigma2=((length(yy)-length(b))/length(yy))*(resid2'*resid2);
varb2=sigma2*inv(xx'*xx);

beta(1,2,1)=b(1,1);  beta(2,2,1)=b(2,1); beta(3,2,1)=b(3,1); beta(4,2,1)=b(4,1); beta(5,2,1)=b(5,1);

p2(:,:,1)=varb2;

Q=sigma2;
R=varb2;

F2=eye(5);
%--------------------------------kalman Filter for Parameter Estimation in the Second Equation of VAR--------------------------------------
for t=1:T-1
    p2(:,:,t+1)=F2*(p2(:,:,t)-(p2(:,:,t)*x(t,:)'*inv(x(t,:)*p2(:,:,t)*x(t,:)'+Q)*x(t,:)*p2(:,:,t)))*F2'+R*R';
    beta(:,2,t+1)=F2*beta(:,2,t)+F2*p2(:,:,t)*x(t,:)'*inv(x(t,:)*p2(:,:,t)*x(t,:)'+Q)*(y(t,2)-x(t,:)*beta(:,2,t));
    betat(:,2,t)=beta(:,2,t)+p2(:,:,t)*x(t,:)'*inv(x(t,:)*p2(:,:,t)*x(t,:)'+Q)*(y(t,2)-x(t,:)*beta(:,2,t));
end

%--------------------------------Parameter Values from OLS in the Third Equation of VAR--------------------------------------

yy=data3(:,3);
xx=lagmatrix([ones(length(yy),1) data3],1);
yy=yy(2:end,:);
xx=xx(2:end,:);
b=inv(xx'*xx)*xx'*yy;
resid3=yy-xx*b;
sigma3=((length(yy)-length(b))/length(yy))*(resid3'*resid3);
varb3=sigma3*inv(xx'*xx);



beta(1,3,1)=b(1,1);   beta(2,3,1)=b(2,1); beta(3,3,1)=b(3,1); beta(4,3,1)=b(4,1); beta(5,3,1)=b(5,1);

p3(:,:,1)=varb3;

Q=sigma3;
R=varb3;

F3=eye(5);
%--------------------------------kalman Filter for Parameter Estimation in the Third Equation of VAR--------------------------------------
for t=1:T-1
    p3(:,:,t+1)=F3*(p3(:,:,t)-(p3(:,:,t)*x(t,:)'*inv(x(t,:)*p3(:,:,t)*x(t,:)'+Q)*x(t,:)*p3(:,:,t)))*F3'+R*R';
    beta(:,3,t+1)=F3*beta(:,3,t)+F3*p3(:,:,t)*x(t,:)'*inv(x(t,:)*p3(:,:,t)*x(t,:)'+Q)*(y(t,3)-x(t,:)*beta(:,3,t));
    betat(:,3,t)=beta(:,3,t)+p3(:,:,t)*x(t,:)'*inv(x(t,:)*p3(:,:,t)*x(t,:)'+Q)*(y(t,3)-x(t,:)*beta(:,3,t));
end


%--------------------------------Parameter Values from OLS in the Fourth Equation of VAR--------------------------------------


yy=data3(:,4);
xx=lagmatrix([ones(length(yy),1) data3],1);
yy=yy(2:end,:);
xx=xx(2:end,:);
b=inv(xx'*xx)*xx'*yy;
resid4=yy-xx*b;
sigma4=((length(yy)-length(b))/length(yy))*(resid4'*resid4);
varb4=sigma4*inv(xx'*xx);


beta(1,4,1)=b(1,1);   beta(2,4,1)=b(2,1);  beta(3,4,1)=b(3,1); beta(4,4,1)=b(4,1); beta(5,4,1)=b(5,1);

p4(:,:,1)=varb4;

Q=sigma4;
R=varb4;

F4=eye(5);
%--------------------------------kalman Filter for Parameter Estimation in the Fourth Equation of VAR--------------------------------------
for t=1:T-1
    p4(:,:,t+1)=F4*(p4(:,:,t)-(p4(:,:,t)*x(t,:)'*inv(x(t,:)*p4(:,:,t)*x(t,:)'+Q)*x(t,:)*p4(:,:,t)))*F4'+R*R';
    beta(:,4,t+1)=F4*beta(:,4,t)+F4*p4(:,:,t)*x(t,:)'*inv(x(t,:)*p4(:,:,t)*x(t,:)'+Q)*(y(t,4)-x(t,:)*beta(:,4,t));
    betat(:,4,t)=beta(:,4,t)+p4(:,:,t)*x(t,:)'*inv(x(t,:)*p4(:,:,t)*x(t,:)'+Q)*(y(t,4)-x(t,:)*beta(:,4,t));
end


%%%%Forecasting 1 to 4 Step Ahead%%%%


% 1-step-ahead forecast
X_tplus1 = [1 y(end,:)];

for w=1:4
    
    onestep(1,w)=X_tplus1*betat(:,w,end);
    
end
if isnan(one_step)~=1
    onestep(1,1)=one_step;
end

% 2-step-ahead forecast

X_tplus2 = [1 onestep(1,:)];



for w=1:4
    
    twostep(1,w)=X_tplus2*betat(:,w,end);
end


%3-step-ahead forecast


X_tplus3 =[1 twostep(1,:)];

for w=1:4
    
    threestep(1,w)=X_tplus3*betat(:,w,end);
end

% 4-step-ahead forecast


X_tplus4 = [1 threestep(1,:)];

for w=1:4
    
    fourstep(1,w)=X_tplus4*betat(:,w,end);
end


forecast=[onestep(1,1) twostep(1,1) threestep(1,1) fourstep(1,1)];




