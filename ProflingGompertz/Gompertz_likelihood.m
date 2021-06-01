function [e] = Gompertz_likelihood(lambda, K, N0, sigma, t,Ndata)
% [t,M]= ode45(@(t,y) lambda*y.*log(K/y), t, N0);
% N=M;
N =  K*exp(log(N0/K)*exp(-1*lambda*t'));
e=0;
y= log(normpdf(N',Ndata,sigma)); 
e=sum(y);