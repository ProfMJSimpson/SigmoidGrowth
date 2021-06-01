clear all
clc
t=[0, 769,1140,1488,1876,2233,2602,2889,3213,3621,4028]; %Observation time
Ndata=[2.352254642,4.396074415,8.434146341,22.25079365,38.9,59.04803013,67.84648814,69.51641791,74.09765494,82.29230769,80.88291457];

beta0  = 1;       %Initial estimate mle(1,1)
lambda0 = 0.0025; %Initial estimate mle(1,2)
K0 = 80.0;        %Initial estimate mle(1,3)
N00=2;            %Initial estimate mle(1,4)
sigma0=100.0;       %Initial estimate mle(1,5)
mle=zeros(1,5); %mle beta, lambda, K, N0
options = optimoptions('fmincon','Display','iter');
nonlcon=[];
gs = GlobalSearch;

%% MLE
funmle=@(n) - Richards_likelihood(n(1), n(2), n(3), n(4),n(5), t, Ndata);
[mle,nLL] = fmincon(funmle,[beta0,lambda0,K0,N00,sigma0],[],[],[],[],[0.0, 0.0,0.0,0.0,0.0],[],nonlcon,options);
% problem = createOptimProblem('fmincon','x0',[beta0,lambda0,K0,N00,sigma0],...
%      'objective',funmle,'lb',[1e-10,1e-10,1e-10,1e-10,1e-10],'ub',[]);
%  x = run(gs,problem);
%Uncomment lines 21-23 to perform the Global search.  Leave commented to
%run the local search only.

Nmle = mle(1,3)*mle(1,4)./(mle(1,4)^mle(1,1)+(mle(1,3)^mle(1,1)-mle(1,4)^mle(1,1)).*exp(-1*mle(1,1)*mle(1,2)*t')).^(1/mle(1,1));
Var = mle(1,5)*ones(size(t));

%Compute R2
R2 = 1.0 - sum((Ndata-[Nmle]').^2)/sum((Ndata-mean(Ndata)).^2);

subplot(2,1,1)
plot(t,Nmle,'g','LineWidth',2)
hold on
plot(t, Ndata,'o','MarkerFaceColor','b','MarkerEdgeColor','b')
xlabel('time')
ylabel('Density')
title('Richards')
legend( sprintf('R^2 %f', R2),'Location','northwest' )


subplot(2,1,2)
res=(Ndata-[Nmle]')/sqrt(mle(1,5));
plot(t,res,'s','MarkerFaceColor','b','MarkerEdgeColor','b')
yline(0,'r','Linewidth',2)
xlabel('time')
ylabel('Normalised Residual')
ylim([-1.2*max(abs(res)) 1.2*max(abs(res))])






