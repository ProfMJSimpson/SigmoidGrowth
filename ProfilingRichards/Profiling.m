clear all
clc
rng(1);
npts=50;
TH=-1.9207;
 t=[0, 769,1140,1488,1876,2233,2602,2889,3213,3621,4028]; %Observation time
 Ndata=[2.352254642,4.396074415,8.434146341,22.25079365,38.9,59.04803013,67.84648814,69.51641791,74.09765494,82.29230769,80.88291457];


beta0=0.3;
lambda0=0.005;
K0=80;
N00=1.0;
sigma0=100.0;




mle=zeros(1,5); %mle: beta, lambda, K, N0, sigma



options = optimoptions('fmincon','Display','iter','MaxIterations',5000);
nonlcon=[];
gs = GlobalSearch;


%% MLE
funmle=@(n) - Richards_Likelihood(n(1), n(2), n(3), n(4), n(5), t, Ndata);
[mle,nLL] = fmincon(funmle,[beta0,lambda0,K0,N00,sigma0],[],[],[],[],[],[],nonlcon,options);
 problem = createOptimProblem('fmincon','x0',[beta0,lambda0,K0,N00,sigma0],...
      'objective',funmle,'lb',[1e-10,1e-10,1e-10,1e-10,1e-10],'ub',[]);
x = run(gs,problem);
%At this point mle is the mle from the local search whereas x is the mle
%from the global search.  

mle=x; %this line accepts the global results but you might want to compare them before proceeding


Nmle= mle(1,3)*mle(1,4)./(mle(1,4)^mle(1,1)+(mle(1,3)^mle(1,1)-mle(1,4)^mle(1,1)).*exp(-1*mle(1,1)*mle(1,2)*t')).^(1/mle(1,1));

figure
plot(t, Ndata,'o')
hold on
plot(t,Nmle,'g')
xlabel('time')
ylabel('Density')



beta_min=0.1; 
beta_max=0.9;
lambda_min=0.001;
lambda_max=0.01;
K_min=70;
K_max=90;
N0_min=0.01;
N0_max=1.0;
sigma_min=0.1;
sigma_max=5.0;
betalambda_min=0.0010;
betalambda_max=0.0030;




%% Code to Profile lambda
rrange = linspace(lambda_min,lambda_max,npts);
rrange=[rrange,mle(1,2)];
rrange=sort(rrange);
nrange=zeros(4,numel(rrange)); %Vector to hold elements of the nuisance parameters
                               %first row is beta, second row is K,
                               %third row is N0, fourth row is eta
                               
lhoodr=zeros(1,numel(rrange));

for i = 1:numel(rrange)
    i;
rr=rrange(i);
funr=@(n) - Richards_Likelihood(n(1), rr, n(2), n(3), n(4), t, Ndata);


      if i==1
      n0=[beta0,K0,N00,sigma0];
      elseif i > 1
      n0=[nrange(:,i-1)];
      end

[nrange(:,i)] = fmincon(funr,n0,[],[],[],[],([0,0,0,0]),[]);
end


for i=1:numel(rrange)
lhoodr(i)=-Richards_Likelihood(nrange(1,i),rrange(1,i),nrange(2,i),nrange(3,i),nrange(4,i),t,Ndata); 
end



lhood2r = min(lhoodr)-lhoodr;
subplot(6,1,1)
plot(rrange,lhood2r,'k','LineWidth',2)
hold on
xline(mle(1,2),'r')
xlim([rrange(1) rrange(numel(rrange))])
ylim([-3, 0])
yline(TH)
xlabel('lambda')
ylabel('Profile')

%Now write a piece of code that will identify the CI
il=1;
ir=numel(rrange)-1;
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
krange = linspace(K_min,K_max,npts);
krange=[krange,mle(1,3)];
krange=sort(krange);
nrange=zeros(4,numel(krange)); %Vector to hold elements of the nuisance parameters
                               %first row is beta, second row is lambda,
                               %third row is N0, fourth row is eta
lhoodk=zeros(1,numel(krange));

for i = 1:numel(krange)
    i
kk=krange(i);
funk=@(n) - Richards_Likelihood(n(1), n(2), kk, n(3), n(4), t, Ndata);


      if i==1
      n0=[beta0,lambda0,N00,sigma0];
      elseif i > 1
      n0=[nrange(:,i-1)];
      end
 
[nrange(:,i)] = fmincon(funk,n0,[],[],[],[],([0,0,0,0]),[]);
end


for i=1:numel(krange)
lhoodk(i)=-Richards_Likelihood(nrange(1,i),nrange(2,i),krange(1,i),nrange(3,i),nrange(4,i),t,Ndata); 
end



lhood2k = min(lhoodk)-lhoodk;
subplot(6,1,2)
plot(krange,lhood2k,'k','LineWidth',2)
hold on
xline(mle(1,3),'r')
xlim([krange(1) krange(numel(krange))])
ylim([-3, 0])
yline(TH)
xlabel('K')
ylabel('Profile')

%Now write a piece of code that will identify the CI
il=1;
ir=numel(krange)-1;
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





%% Code to Profile beta
brange = linspace(beta_min,beta_max,npts);
brange=[brange,mle(1,1)];
brange=sort(brange);
nrange=zeros(4,numel(brange)); %Vector to hold elements of the nuisance parameters
                               %first row is lambda, second row is K,
                               %third row is N0, fourth row is eta0

lhoodb=zeros(1,numel(brange));

for i =numel(brange):-1:1
    i
b=brange(i);
funb=@(n) - Richards_Likelihood(b, n(1), n(2), n(3), n(4), t, Ndata);


      if i==numel(brange)
      n0=[lambda0, K0, N00, sigma0];
      elseif i > 1
      n0=[nrange(:,i+1)];
      end

%[nrange(:,i)] = fmincon(funb,n0,[],[],[],[],([1e-10,1e-10,1e-10,1e-10]),[]);
problem = createOptimProblem('fmincon','x0',[lambda0,K0,N00,sigma0],...
      'objective',funb,'lb',[1e-10,1e-10,1e-10,1e-10],'ub',[]);
[nrange(:,i)]= run(gs,problem);
 
 
end


for i=1:numel(brange)
lhoodb(i)=-Richards_Likelihood(brange(1,i),nrange(1,i),nrange(2,i),nrange(3,i),nrange(4,i),t,Ndata); 
end



lhood2b = min(lhoodb)-lhoodb;
subplot(6,1,3)
plot(brange,lhood2b,'k','LineWidth',2)
hold on
xline(mle(1,1),'r')
xlim([brange(1) brange(numel(brange))])
ylim([-3, 0])
yline(TH)
xlabel('beta')
ylabel('Profile')

%Now write a piece of code that will identify the CI
il=1;
ir=numel(brange)-1;
for i=1:numel(brange)-1
    if lhood2b(i) < TH && lhood2b(i+1) > TH
        il=i;
    elseif lhood2b(i) > TH && lhood2b(i+1) < TH
        ir=i;
    end
end

%Approximate the CI
CI_lower_b=(TH*(brange(il+1)-brange(il))+lhood2b(il+1)*brange(il)-lhood2b(il)*brange(il+1))/(lhood2b(il+1)-lhood2b(il)); %this expression linearly interpolates
CI_upper_b=(TH*(brange(ir+1)-brange(ir))+lhood2b(ir+1)*brange(ir)-lhood2b(ir)*brange(ir+1))/(lhood2b(ir+1)-lhood2b(ir)); %this expression linearly interpolates
xline(CI_lower_b,'g')
xline(CI_upper_b,'c')






%% Code to Profile N0
N0range = linspace(N0_min,N0_max,npts);
N0range=[N0range,mle(1,4)];
N0range=sort(N0range);
nrange=zeros(4,numel(N0range)); %Vector to hold elements of the nuisance parameters
                               %first row is beta, second row is lambda,
                               %third row is K, fourth row is eta
lhoodN0=zeros(1,numel(N0range));

for i = 1:numel(N0range)
N0=N0range(i);
i
funN0=@(n) - Richards_Likelihood(n(1), n(2), n(3), N0, n(4), t, Ndata);

      if i==1
      n0=[beta0,lambda0,K0,sigma0];
      elseif i > 1
      n0=[nrange(:,i-1)];
      end
 
[nrange(:,i)] = fmincon(funN0,n0,[],[],[],[],([0,0,0,0]),[]);
end


for i=1:numel(N0range)
lhoodN0(i)=-Richards_Likelihood(nrange(1,i),nrange(2,i),nrange(3,i),N0range(1,i),nrange(4,i),t,Ndata); 
end


lhood2N0 = min(lhoodN0)-lhoodN0;
subplot(6,1,4)
plot(N0range,lhood2N0,'k','LineWidth',2)
hold on
xline(mle(1,4),'r')
xlim([N0range(1) N0range(numel(N0range))])
ylim([-3, 0])
yline(TH)
xlabel('N(0)')
ylabel('Profile')

%Now write a piece of code that will identify the CI
il=1;
ir=numel(N0range)-1;
for i=1:numel(N0range)-1
    if lhood2N0(i) < TH && lhood2N0(i+1) > TH
        il=i;
    elseif lhood2N0(i) > TH && lhood2N0(i+1) < TH
        ir=i;
    end
end

%Approximate the CI
if il > 1 && ir < 100
CI_lower_N0=(TH*(N0range(il+1)-N0range(il))+lhood2N0(il+1)*N0range(il)-lhood2N0(il)*N0range(il+1))/(lhood2N0(il+1)-lhood2N0(il)); %this expression linearly interpolates
CI_upper_N0=(TH*(N0range(ir+1)-N0range(ir))+lhood2N0(ir+1)*N0range(ir)-lhood2N0(ir)*N0range(ir+1))/(lhood2N0(ir+1)-lhood2N0(ir)); %this expression linearly interpolates
else
CI_lower_N0=N0_min;
CI_upper_N0=N0_max;
end

xline(CI_lower_N0,'g')
xline(CI_upper_N0,'c')


%% Code to Profile sigma
sigmarange=linspace(sigma_min,sigma_max,npts);
sigmarange=[sigmarange,mle(1,5)];
sigmarange=sort([sigmarange]);
nrange=zeros(4,numel(sigmarange)); %Vector to hold elements of the nuisance parameters
                               %first row is beta, second row is lambda,
                               %third row is K, N0
lhoodsigma=zeros(1,numel(sigmarange));

for i = numel(sigmarange):-1:1
sg=sigmarange(i);
i
funsigma=@(n) - Richards_Likelihood(n(1), n(2), n(3), n(4), sg, t, Ndata);


      if i==numel(sigmarange)
      n0=[mle(1,1), mle(1,2), mle(1,3), mle(1,4)];
      %n0=[lambda0, K0, N00];
      elseif i > 1
      n0=[nrange(:,i+1)];
      end
 
[nrange(:,i)] = fmincon(funsigma,n0,[],[],[],[],([0,0,0,0]),[]);
end


for i=1:numel(sigmarange)
lhoodsigma(i)=-Richards_Likelihood(nrange(1,i),nrange(2,i),nrange(3,i),nrange(4,i),sigmarange(1,i),t,Ndata); 
end



lhood2sigma = min(lhoodsigma)-lhoodsigma;
subplot(6,1,5)
plot(sigmarange,lhood2sigma,'k','LineWidth',2)
hold on
xline(mle(1,5),'r')
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






%% Code to profile product of beta*lambda
brrange = linspace(betalambda_min,betalambda_max,npts);
brrange=[brrange,mle(1,1)*mle(1,2)];
brrange=sort(brrange);

%Vector to hold elements of the nuisance parameters
%(lambda,K,N0,sigma)
                               
lhoodr=zeros(1,numel(brrange));

for i = 1:numel(brrange)
    i
brr=brrange(i);
funbr=@(n) - Richards_Likelihood_reparmeterised(brr, n(1), n(2), n(3), n(4), t, Ndata);


      if i==1
      n0=[lambda0,K0,N00,sigma0];
      elseif i > 1
      n0=[nrange(:,i-1)];
      end

%[nrange(:,i)] = fmincon(funbr,n0,[],[],[],[],([0,0,0,0]),[]);
problem = createOptimProblem('fmincon','x0',[lambda0,K0,N00,sigma0],...
      'objective',funbr,'lb',[1e-10,1e-10,1e-10,1e-10],'ub',[]);
[nrange(:,i)]= run(gs,problem);

end


for i=1:numel(brrange)
lhoodbr(i)=-Richards_Likelihood_reparmeterised(brrange(1,i),nrange(1,i),nrange(2,i),nrange(3,i),nrange(4,i),t,Ndata); 
end



lhood2br = min(lhoodbr)-lhoodbr;
subplot(6,1,6)
plot(brrange,lhood2br,'k','LineWidth',2)
xlim([brrange(1) brrange(numel(brrange))])
ylim([-3, 0])
yline(TH)
xlabel('beta times lambda')
ylabel('Profile')

%Now write a piece of code that will identify the CI
il=1;
ir=numel(brrange)-1;
for i=1:numel(brrange)-1
    if lhood2br(i) < TH && lhood2br(i+1) > TH
        il=i;
    elseif lhood2br(i) > TH && lhood2br(i+1) < TH
        ir=i;
    end
end


%Approximate the CI
CI_lower_br=(TH*(brrange(il+1)-brrange(il))+lhood2br(il+1)*brrange(il)-lhood2br(il)*brrange(il+1))/(lhood2br(il+1)-lhood2br(il)); %this expression linearly interpolates
CI_upper_br=(TH*(brrange(ir+1)-brrange(ir))+lhood2br(ir+1)*brrange(ir)-lhood2br(ir)*brrange(ir+1))/(lhood2br(ir+1)-lhood2br(ir)); %this expression linearly interpolates
xline(CI_lower_br,'g')
xline(CI_upper_br,'c')


 %% Code to construct bivariate profile (beta,lambda)
  %mle: beta, lambda, K, N0, sigma
  
 lambda_min=0.001;
 lambda_max=0.010;
 beta_min=0.1;
 beta_max=1.0;
  
 rrange = linspace(lambda_min,lambda_max,npts);
 brange=linspace(beta_min,beta_max,npts);
 
 nrange=zeros(3,numel(rrange),numel(brange)); %Vector to hold elements of the two nuisance parameters
                                              %first row is K,
                                              %second row is N0,
                                              %third row is sigma
 
 lhoodrb=zeros(numel(rrange),numel(brange));
 
 for i = 1:numel(rrange)
     for j=1:numel(brange)
 r=rrange(i);
 b=brange(j);
 i
 funrb=@(n) -Richards_Likelihood(b, r, n(1), n(2), n(3), t, Ndata);
 
        
       n0=[K0,N00,sigma0];
       
 nrange(:,i,j) = fmincon(funrb,n0,[],[],[],[],[0,0,0],[]);
     end
 end
 
 
 for i=1:numel(rrange)
    for j=1:numel(brange) 
 lhoodrb(i,j)=-Richards_Likelihood(brange(1,j),rrange(1,i),nrange(1,i,j),nrange(2,i,j),nrange(3,i,j),t,Ndata);
    end
 end
 
 
 
lhood2rb = nLL-lhoodrb;
 figure
 
 contourf(rrange,brange,lhood2rb,'LevelList',[-3.00],'Linewidth',2);
 
 %colorbar
 hold on
 plot(mle(1,2),mle(1,1),'ro','MarkerFaceColor','r')
 xlabel("lambda");
 ylabel("beta");
 
 
 
 
