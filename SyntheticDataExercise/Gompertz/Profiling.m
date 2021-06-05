clear all
clc
rng(1);
bb= [0.05]; %linspace(0.05,1.5,30);
npts=100;
TH=-1.9207;
err =zeros(1,numel(bb));
width=zeros(1, numel(bb));


for ii=1:length(bb)
beta =  bb(ii);
lambda = 0.0055;
K=80;
N0=1.0;
NMes=20;
omega = 0.999; % This is the proportion of time we examine
sigma=2.0/1.;
tt = log((K^beta-N0^beta)/((N0/omega)^beta-N0^beta))/(lambda*beta);
t = linspace(0,tt,NMes);
Ndata = K*N0./(N0^beta+(K^beta-N0^beta).*exp(-1*lambda*beta*t)).^(1/beta);
subplot(2,1,1)
  plot(t*beta*lambda,Ndata/K,'r','LineWidth',1)
Ndata=Ndata+normrnd(0,sigma,[1,length(t)]);
 hold on
 plot(t*beta*lambda,Ndata/K,'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',2)
%  ylim([0 1.2])



K0=K;
N00=N0;
lambda0=lambda;
sigma0=100.0;






mle=zeros(1,4); %mle lambda, K, N0, eta

options = optimoptions('fmincon','Display','iter','MaxIterations',5000);
nonlcon=[];
gs = GlobalSearch;


%% MLE
funmle=@(n) - Gompertz_likelihood(n(1), n(2), n(3), n(4), t, Ndata);
%[mle,nLL] = fmincon(funmle,[lambda0,K0,N00,sigma0],[],[],[],[],[0,0,0,0],[],nonlcon,options);
problem = createOptimProblem('fmincon','x0',[lambda0,K0,N00,sigma0],'objective',funmle,'lb',[1e-10,1e-10,1e-10,1e-10],'ub',[]);
x = run(gs,problem);
mle=x; %this line accepts the global results but you might want to compare them before proceeding


err(ii) = beta*lambda - mle(1,1);




Nmle=mle(1,2)*exp(log(mle(1,3)/mle(1,2))*exp(-1*mle(1,1)*t'));

hold on
plot(beta*t*lambda,Nmle/mle(1,2),'g')
xlabel('time')
ylabel('Density')


%Calculate and plot the scaled-residuals
Var = mle(1,4)*ones(size(t));



res=(Ndata-[Nmle]')/sqrt(mle(1,4));
subplot(2,1,2)
plot(t*lambda*beta,res,'o','MarkerFaceColor','g','MarkerEdgeColor','g')
yline(0,'r','Linewidth',2)
xlabel('time')
ylabel('Normalised Residual')
ymax = max(abs(res));
ylim([-1.2*ymax 1.2*ymax])




 lambda_min=mle(1,1)/2; 
 lambda_max=mle(1,1)*2;

 
%% Code to Profile lambda
rrange = linspace(lambda_min,lambda_max,npts);
rrange=[rrange,mle(1,1)];
rrange=sort(rrange);

nrange=zeros(3,numel(rrange)); %Vector to hold elements of thenuisance parameters
                               %first row is K, second row is N0,
                               %third row is eta
lhoodr=zeros(1,numel(rrange));

for i = 1:numel(rrange)
rr=rrange(i);
funr=@(n) - Gompertz_likelihood(rr, n(1), n(2), n(3), t, Ndata);


      if i==1
      n0=[K0, N00, sigma0];
      elseif i > 1
      n0=[nrange(:,i-1)];
      end

%[nrange(:,i)] = fmincon(funr,n0,[],[],[],[],([0,0,0]),[]);
problem = createOptimProblem('fmincon','x0',[n0],'objective',funr,'lb',[1e-10,1e-10,1e-10],'ub',[]);
[nrange(:,i)] = run(gs,problem);
end


for i=1:numel(rrange)
lhoodr(i)=- Gompertz_likelihood(rrange(1,i),nrange(1,i),nrange(2,i),nrange(3,i),t,Ndata); 
end


 figure
lhood2r = min(lhoodr)-lhoodr;
plot(rrange,lhood2r,'g','LineWidth',2)
hold on
xline(mle(1,1),'r')
xline(beta*lambda,'c')
ylim([-3, 0])
xlim([0 lambda_max])
yline(TH)
xlabel('lambda')
ylabel('Profile')

%Now write a piece of code that will identify the CI
il=1;
ir=numel(rrange);
for i=1:numel(rrange)-1
    if lhood2r(i) < TH && lhood2r(i+1) > TH
        il=i;
    elseif lhood2r(i) > TH && lhood2r(i+1) < TH
        ir=i;
    end
end

%Approximate the lower CI
CI_lower_lambda=(TH*(rrange(il+1)-rrange(il))+lhood2r(il+1)*rrange(il)-lhood2r(il)*rrange(il+1))/(lhood2r(il+1)-lhood2r(il)); %this expression linearly interpolates
CI_upper_lambda=(TH*(rrange(ir+1)-rrange(ir))+lhood2r(ir+1)*rrange(ir)-lhood2r(ir)*rrange(ir+1))/(lhood2r(ir+1)-lhood2r(ir)); %this expression linearly interpolates

% xline(CI_lower_lambda,'g')
% xline(CI_upper_lambda,'c')
CI_width_lambda = CI_upper_lambda-CI_lower_lambda;
width(1,ii) = CI_width_lambda;


end

figure
subplot(1,2,1)
plot(bb,err,'og','MarkerFaceColor','g','MarkerEdgeColor','g')
xlabel('\beta')
ylabel('\beta \lambda_3 - \lambda_1')
subplot(1,2,2)
plot(bb,width,'og','MarkerFaceColor','g','MarkerEdgeColor','g')
xlabel('\beta')
ylabel('CI Width')
