function [e] = Richards_Gaussianlikelihood(beta, lambda, K, N0, sigma, t,Ndata)
N = K*N0./(N0^beta+(K^beta-N0^beta).*exp(-1*beta*lambda*t')).^(1/beta);
e=0;
y= log(normpdf(N',Ndata,sigma)); 
e=sum(y);