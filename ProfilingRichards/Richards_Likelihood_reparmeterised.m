function [e] = Richards_Likelihood_reparmeterised(betalambda, lambda, K, N0, sigma, t,Ndata)
beta = betalambda/lambda;
N = K*N0./(N0^beta+(K^beta-N0^beta).*exp(-1*beta*lambda*t')).^(1/beta);
e=0;
y= log(normpdf(N',Ndata,sigma)); 
e=sum(y);