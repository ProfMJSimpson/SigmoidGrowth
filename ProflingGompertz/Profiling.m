clear all
clc
rng(1);
npts=100;
TH=-1.9207;
t=[0, 769,1140,1488,1876,2233,2602,2889,3213,3621,4028]; %Observation time
Ndata=[2.352254642,4.396074415,8.434146341,22.25079365,38.9,59.04803013,67.84648814,69.51641791,74.09765494,82.29230769,80.88291457];


lambda0=0.002;
K0=80;
N00=1.0;
sigma0=10.0;






mle=zeros(1,4); %mle lambda, K, N0, sigma

options = optimoptions('fmincon','Display','iter','MaxIterations',5000);
nonlcon=[];
gs = GlobalSearch;


%% MLE
funmle=@(n) - Gompertz_likelihood(n(1), n(2), n(3), n(4), t, Ndata);
[mle,nLL] = fmincon(funmle,[lambda0,K0,N00,sigma0],[],[],[],[],[0,0,0,0],[],nonlcon,options);
% problem = createOptimProblem('fmincon','x0',[lambda0,K0,N00,sigma0],'objective',funmle,'lb',[1e-10,1e-10,1e-10,1e-10],'ub',[]);
% x = run(gs,problem);
%mle=x; 

Nmle=mle(1,2)*exp(log(mle(1,3)/mle(1,2))*exp(-1*mle(1,1)*t'));
figure
plot(t, Ndata,'o')
hold on
plot(t,Nmle,'g')
xlabel('time')
ylabel('Density')



lambda_min=0.0010;
lambda_max=0.0050;
K_min=70.0;
K_max=90;
N0_min=0.000001;
N0_max=0.015;
sigma_min=1.0;
sigma_max=5.0;







%% Code to Profile lambda
rrange = linspace(lambda_min,lambda_max,npts);
rrange=[rrange,mle(1,1)];
rrange=sort(rrange);

nrange=zeros(3,numel(rrange)); %Vector to hold elements of thenuisance parameters
                               %first row is K, second row is N0,
                               %third row is eta
lhoodr=zeros(1,numel(rrange));

for i = numel(rrange):-1:1
    i
rr=rrange(i);
funr=@(n) - Gompertz_likelihood(rr, n(1), n(2), n(3), t, Ndata);


      if i==numel(rrange)
      n0=[K0, N00, sigma0];
      %n0=[mle(1,2), mle(1,3), mle(1,4)]
      elseif i > 1
      n0=[nrange(:,i+1)];
      end

%[nrange(:,i)] = fmincon(funr,n0,[],[],[],[],([0,0,0]),[]);
problem = createOptimProblem('fmincon','x0',[n0],'objective',funr,'lb',[1e-10,0,0],'ub',[]);
[nrange(:,i)] = run(gs,problem);
end


for i=1:numel(rrange)
lhoodr(i)=- Gompertz_likelihood(rrange(1,i),nrange(1,i),nrange(2,i),nrange(3,i),t,Ndata); 
end



lhood2r = min(lhoodr)-lhoodr;
subplot(4,1,1)
plot(rrange,lhood2r,'k','LineWidth',2)
hold on
xline(mle(1,1),'r')
xlim([rrange(1) rrange(numel(rrange))])
ylim([-3, 0])
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

%Approximate the CI
CI_lower_lambda=(TH*(rrange(il+1)-rrange(il))+lhood2r(il+1)*rrange(il)-lhood2r(il)*rrange(il+1))/(lhood2r(il+1)-lhood2r(il)); %this expression linearly interpolates
CI_upper_lambda=(TH*(rrange(ir+1)-rrange(ir))+lhood2r(ir+1)*rrange(ir)-lhood2r(ir)*rrange(ir+1))/(lhood2r(ir+1)-lhood2r(ir)); %this expression linearly interpolates
xline(CI_lower_lambda,'g')
xline(CI_upper_lambda,'c')



%% Code to Profile K
krange = linspace(70,90,npts);
krange=[krange,mle(1,2)];
krange=sort(krange);
nrange=zeros(3,numel(krange)); %Vector to hold elements of the nuisance parameters
                               %first row is lambda, second row is N0,
                               %third row is eta0


for i = numel(krange):-1:1
    i
kk=krange(i);
funk=@(n) -Gompertz_likelihood(n(1), kk, n(2), n(3), t, Ndata);


      if i==numel(krange)
      n0=[lambda0, N00, sigma0];
      %n0=[mle(1,1),mle(1,3),mle(1,4)];
      elseif i > 1
      n0=[nrange(:,i+1)];
      end
 
[nrange(:,i)] = fmincon(funk,n0,[],[],[],[],([0,0,0]),[]);
%problem = createOptimProblem('fmincon','x0',[n0],'objective',funk,'lb',[1e-10,1e-10,1e-10],'ub',[]);
%[nrange(:,i)] = run(gs,problem);

end


for i=1:numel(krange)
lhoodk(i)=- Gompertz_likelihood(nrange(1,i),krange(1,i),nrange(2,i),nrange(3,i),t,Ndata); 
end

lhood2k = min(lhoodk)-lhoodk;
subplot(4,1,2)
plot(krange,lhood2k,'k','LineWidth',2)
hold on
xline(mle(1,2),'r')
xlim([krange(1) krange(numel(krange))])
ylim([-3, 0])
yline(TH)
xlabel('K')
ylabel('Profile')



%Now write a piece of code that will identify the CI
il=1;
ir=numel(krange);
for i=1:numel(krange)-1
    if lhood2k(i) < TH && lhood2k(i+1) > TH
        il=i;
    elseif lhood2k(i) > TH && lhood2k(i+1) < TH
        ir=i;
    end
end

%Approximate the CI
CI_lower_K=(TH*(krange(il+1)-krange(il))+lhood2k(il+1)*krange(il)-lhood2k(il)*krange(il+1))/(lhood2k(il+1)-lhood2k(il)); %this expression linearly interpolates
CI_upper_K=(TH*(krange(ir+1)-krange(ir))+lhood2k(ir+1)*krange(ir)-lhood2k(ir)*krange(ir+1))/(lhood2k(ir+1)-lhood2k(ir)); %this expression linearly interpolates
xline(CI_lower_K,'g')
xline(CI_upper_K,'c')


%% Code to Profile N0
N0range = linspace(N0_min,N0_max,npts);
N0range=[N0range,mle(1,3)];
N0range=sort(N0range);
nrange=zeros(3,numel(N0range)); %Vector to hold elements of the nuisance parameters
                               %first row is lambda, second row is K,
                               %third row is eta
lhoodN0=zeros(1,numel(N0range));

for i = 1:numel(N0range)
N0=N0range(i);
i
funN0=@(n) - Gompertz_likelihood(n(1), n(2), N0, n(3), t, Ndata);


      if i==1
      n0=[lambda0, K0, sigma0];
      elseif i > 1
      n0=[nrange(:,i-1)];
      end
 
[nrange(:,i)] = fmincon(funN0,n0,[],[],[],[],([0,0,0]),[]);
%problem = createOptimProblem('fmincon','x0',[n0],'objective',funN0,'lb',[1e-10,1e-10,1e-10],'ub',[]);
%[nrange(:,i)] = run(gs,problem);
end


for i=1:numel(N0range)
lhoodN0(i)=-Gompertz_likelihood(nrange(1,i),nrange(2,i),N0range(1,i),nrange(3,i),t,Ndata); 
end



lhood2N0 = min(lhoodN0)-lhoodN0;
subplot(4,1,3)
plot(N0range,lhood2N0,'k','LineWidth',2)
hold on
xline(mle(1,3),'r')
xlim([N0range(1) N0range(numel(N0range))])
ylim([-3, 0])
yline(TH)
xlabel('N(0)')
ylabel('Profile')



%Now write a piece of code that will identify the CI
il=1;
ir=numel(N0range);
for i=1:numel(N0range)-1
    if lhood2N0(i) < TH && lhood2N0(i+1) > TH
        il=i;
    elseif lhood2N0(i) > TH && lhood2N0(i+1) < TH
        ir=i;
    end
end

%Approximate the CI
CI_lower_N0=(TH*(N0range(il+1)-N0range(il))+lhood2N0(il+1)*N0range(il)-lhood2N0(il)*N0range(il+1))/(lhood2N0(il+1)-lhood2N0(il)); %this expression linearly interpolates
CI_upper_N0=(TH*(N0range(ir+1)-N0range(ir))+lhood2N0(ir+1)*N0range(ir)-lhood2N0(ir)*N0range(ir+1))/(lhood2N0(ir+1)-lhood2N0(ir)); %this expression linearly interpolates
xline(CI_lower_N0,'g')
xline(CI_upper_N0,'c')


%% Code to Profile sigma
sigmarange=linspace(sigma_min,sigma_max,npts);
sigmarange=[sigmarange,mle(1,4)];
sigmarange=sort([sigmarange]);
nrange=zeros(3,numel(sigmarange)); %Vector to hold elements of the nuisance parameters
                               %first row is lambda, second row is K,
                               %third row is N0
lhoodsigma=zeros(1,numel(sigmarange));

for i = numel(sigmarange):-1:1
sg=sigmarange(i);
i
funsigma=@(n) - Gompertz_likelihood(n(1), n(2), n(3), sg, t, Ndata);


      if i==numel(sigmarange)
      n0=[mle(1,1), mle(1,2), mle(1,3)];
      %n0=[lambda0, K0, N00];
      elseif i > 1
      n0=[nrange(:,i+1)];
      end
 
[nrange(:,i)] = fmincon(funsigma,n0,[],[],[],[],([0,0,0]),[]);
%problem = createOptimProblem('fmincon','x0',[lambda0,K0,sigma0],'objective',funN0,'lb',[1e-10,1e-10,1e-10],'ub',[]);
%[nrange(:,i)] = run(gs,problem);
end


for i=1:numel(sigmarange)
lhoodsigma(i)=-Gompertz_likelihood(nrange(1,i),nrange(2,i),nrange(3,i),sigmarange(1,i),t,Ndata); 
end



lhood2sigma = min(lhoodsigma)-lhoodsigma;
subplot(4,1,4)
plot(sigmarange,lhood2sigma,'k','LineWidth',2)
hold on
xline(mle(1,4),'r')
xlim([sigmarange(1) sigmarange(numel(N0range))])
ylim([-3, 0])
yline(-1.9207)
xlabel('sigma')
ylabel('Profile')


%Now write a piece of code that will identify the CI
for i=1:numel(sigmarange)-1
    if lhood2sigma(i) < TH && lhood2sigma(i+1) > TH
        il=i;
    elseif lhood2sigma(i) > TH && lhood2sigma(i+1) < TH
        ir=i;
    end
end

%Approximate the CI
CI_lower_sigma=(TH*(sigmarange(il+1)-sigmarange(il))+lhood2sigma(il+1)*sigmarange(il)-lhood2sigma(il)*sigmarange(il+1))/(lhood2sigma(il+1)-lhood2sigma(il)); %this expression linearly interpolates
CI_upper_sigma=(TH*(sigmarange(ir+1)-sigmarange(ir))+lhood2sigma(ir+1)*sigmarange(ir)-lhood2sigma(ir)*sigmarange(ir+1))/(lhood2sigma(ir+1)-lhood2sigma(ir)); %this expression linearly interpolates
xline(CI_lower_sigma,'g')
xline(CI_upper_sigma,'c')


