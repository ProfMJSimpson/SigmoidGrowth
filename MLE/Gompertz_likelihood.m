function [e] = Gompertz_likelihood(lambda, K, N0, sigma, t,Ndata)
N =  K*exp(log(N0/K)*exp(-1*lambda*t'));
e=0;
y= log(normpdf(N',Ndata,sigma)); 
e=sum(y);