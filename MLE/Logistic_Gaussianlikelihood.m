function [e] = Logistic_Gaussianlikelihood(lambda, K, N0, sigma, t, Ndata)
N = K*N0./(N0+(K-N0).*exp(-1*lambda*t'));
e=0;
y= log(normpdf(N',Ndata,sigma)); 
e=sum(y);