clear all
clc
rng(1);
%Synthetic data
beta = 1/10.0;
lambda3= 0.01;
K=100;
N0=1.0;
NMes=100;
omega = 0.9999; % This is the proportion of time we examine
tt = log((K^beta-N0^beta)/((N0/omega)^beta-N0^beta))/(lambda3*beta);
t = linspace(0,tt,NMes);
NdataR = K*N0./(N0^beta+(K^beta-N0^beta).*exp(-1*lambda3*beta*t)).^(1/beta);
subplot(1,4,1)
hold on
plot(beta*lambda3*t,NdataR/K,'r','LineWidth',2)




lambda2=lambda3*beta;
NdataG=K.*exp(log(N0/K).*exp(-lambda2.*t));
hold on
plot(lambda3*beta*t,NdataG/K,'g--','LineWidth',2)
%legend('Richards','Gompertz')

lambda1=lambda3*beta;
NdataL=K*N0./(N0+(K-N0).*exp(-lambda1.*t));
plot(lambda3*beta*t,NdataL/K,'b--','LineWidth',2)

xlim([0 10])
ylim([0 1.2])
xticks([0 5 10])
yticks([0 0.5 1.0])
axis square


