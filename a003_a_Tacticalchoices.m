%% 3 Monthly tactical choices
clear
clc
load data.mat

w=60;   % 5 years

%% Group creation
% groups are created considering the period before the current analysis
% window (returns before 2013)
figure
B = corr(returns(1:size(returns,1)-w,:));
Y = pdist(B);
Ylink = linkage(Y);
[h, nodes] = dendrogram(Ylink,'Labels',labels);
clear B Y

% let's see if the group clustering changes in periods of overall negative
% returns. During crisis assets correlation tends to increase
figure
B = corr(returns(benchReturns(1:size(returns,1)-w)<0,:));
Y = pdist(B);
Ylink = linkage(Y);
dendrogram(Ylink,'Labels',labels, 'ColorThreshold', 1.12);
cophenet(Ylink,Y)
clear B Y

% Four groups:
% 1) CZ, HN, PO      2) ES, IT, IR, PT
% 3) UK, NW, SW, SD  4) the rest
PIIS = zeros(18,1)'; PIIS(labels == "PT" | labels == "IT" | labels == "ES" | labels == "IR") = 1;
EAST = zeros(18,1)'; EAST(labels == "CZ" | labels == "PO" | labels == "HN") = 1;
NEU = zeros(18,1)'; NEU(labels == "UK" | labels == "NW" | labels == "SD" | labels == "SW") = 1;
CORE = ones(18,1)' - EAST - PIIS - NEU;

%% Restrict quantities
% Analysis will be restricted to the last 5 years
% Restrict returns to the last 5 years plus initial estimation window w
analysisWindow = 60; % 5 years
returns = returns(end-analysisWindow-w+1:end,:);
benchReturns = benchReturns(end-analysisWindow-w+1:end,:);
% Restrict dates to last 5 years of data
dates = dates(end-analysisWindow+1:end);

% Detailed date range for bigger plots
fineDateRange = (dates(1) : 180 : dates(end));

[r,c]=size(returns);

%% CASE 1 no short selling

% estimation of expected returns
ErS=zeros(r-w,c);   % sample estimator
ErE=zeros(r-w,c);   % equilibrium return
% estimation of covariances
EvS=zeros(r-w,c,c);     % sample estimator
EvE=zeros(r-w,c,c);     % equilibrium
EvW=zeros(r-w,c,c);     % EWMA
l=0.95;                 % smoothing factor for EWMA
EvW0=cov(returns(r-1-w:r-2,:));% initialization for EWMA

for j=1:r-w
    % RETURNS
    % sample estimator
    ErS(j,:)=mean(returns(j:j+w-1,:));
    % equilibrium return - first estimate coefficients and then compute
    % expected returns and allocate
    Y=returns(j:j+w-1,:);
    X=[ones(w,1) benchReturns(j:j+w-1,:)];
    B=(X'*X)\(X'*Y);
    ErE(j,:)=B(2,:)*mean(benchReturns(j:j+w-1,:));
    
    % COVARIANCES
    % sample estimator
    EvS(j,:,:)=cov(returns(j:j+w-1,:));
    % EWMA
    EvW(j,:,:)=l*squeeze(EvW0)+(1-l)*(returns(j+w-1,:)'*returns(j+w-1,:));
    EvW0=EvW(j,:,:);
    % equilibrium estimator
    res = Y - B(2,:) .* benchReturns(j:j+w-1,:);
    EvE(j,:,:) = (B(2,:)')*B(2,:)*var(benchReturns(j:j+w-1,:))+diag(diag(cov(res)));
end
clear Y X B;

% check that matrices are positive definite
cSample = 0; 
cEquilibrium = 0;
cEWMA = 0;
for j=1:size(EvS,1)
    cSample = all(eig(squeeze(EvS(j,:,:)))>eps) + cSample;
    cEquilibrium = all(eig(squeeze(EvE(j,:,:)))>eps) + cEquilibrium;
    cEWMA = all(eig(squeeze(EvW(j,:,:)))>eps) + cEWMA;
end
if (cSample ~= analysisWindow || cEquilibrium ~= analysisWindow || cEWMA ~= analysisWindow)
    disp('Some varcov matrix is not positive definite')
end
% all matrices are positive definite
clear cSample cEquilibrium cEWMA

% evaluation of portfolio weights and realized returns
% loop for evaluation of realized returns in the EW case
% return computation
PortRetEW=sum(returns(w+1:r,:),2)./c;
CREW=cumprod(PortRetEW(:,1)/100+1);


% SAMPLE RETURNS AND COVARIANCES
PortRetS1=zeros(r-w,2); %returns
PortWS1=zeros(2,r-w,c); %weights
% creating no shorting portfolio objects
p=Portfolio;
p=p.setAssetList(labels);
p=p.setDefaultConstraints;
pwgtMS=zeros(c,1);

% EQUILIBRIUM RETURNS AND EWMA COVARIANCES
PortRetE1=zeros(r-w,2); %returns
PortWE1=zeros(2,r-w,c); %weights
p1=Portfolio;                       
p1=p1.setAssetList(labels);
p1=p1.setDefaultConstraints;
pwgtMS1=zeros(c,1);

% EQUILIBRIUM RETURNS AND COVARIANCES
PortRetQ1=zeros(r-w,2); %returns
PortWQ1=zeros(2,r-w,c); %weights
p2=Portfolio;                       
p2=p2.setAssetList(labels);
p2=p2.setDefaultConstraints;
pwgtMS2=zeros(c,1);
for j=1:r-w
    % SAMPLE RETURNS AND COVARIANCES
    % set moments
    p=p.setAssetMoments(ErS(j,:),squeeze(EvS(j,:,:)));
    % GMV
    [PortWS1(1,j,:), PortRetS1(j,1)]=funGMV(p,returns,j,w);
    % MS
    [PortWS1(2,j,:), PortRetS1(j,2), pwgtMS]=funMS(p,returns,j,w,pwgtMS,ErS);
    
    % EQUILIBRIUM RETURNS AND EWMA COVARIANCES
    % set moments
    p1=p1.setAssetMoments(ErE(j,:),squeeze(EvW(j,:,:)));
    % GMV
    [PortWE1(1,j,:), PortRetE1(j,1)]=funGMV(p1,returns,j,w);
    % MS
    [PortWE1(2,j,:), PortRetE1(j,2),pwgtMS1]=funMS(p1,returns,j,w,pwgtMS1,ErE);
    
    % EQUILIBRIUM RETURNS AND COVARIANCES
    % set moments
    p2=p2.setAssetMoments(ErE(j,:),squeeze(EvE(j,:,:)));
    % GMV
    [PortWQ1(1,j,:), PortRetQ1(j,1)]=funGMV(p2,returns,j,w);
    % MS
    [PortWQ1(2,j,:), PortRetQ1(j,2),pwgtMS1]=funMS(p2,returns,j,w,pwgtMS2,ErE);
end

%% area plots for weights over time
figure;
subplot(2,3,1);
area(dates, squeeze(PortWS1(1,:,:)));
% For finer control of dates range:
% set(gca, 'XTick', (dates(1) : 180 : dates(end)) );
title('GMV with sample moments')
ylim([0 1]);
subplot(2,3,4);
area(dates, squeeze(PortWS1(2,:,:)));
title('MS with sample moments')
ylim([0 1]);
subplot(2,3,2);
area(dates, squeeze(PortWE1(1,:,:)));
title('GMV with EWMA')
ylim([0 1]);
subplot(2,3,5);
area(dates, squeeze(PortWE1(2,:,:)));
title('MS with EWMA')
ylim([0 1]);
subplot(2,3,3);
area(dates, squeeze(PortWQ1(1,:,:)));
title('GMV with EQ')
ylim([0 1]);
subplot(2,3,6);
area(dates, squeeze(PortWQ1(2,:,:)));
title('MS with EQ')
ylim([0 1]);
legend(labels)
sgtitle("Weights | No-shorting");

% aggregated group weights
figure;
subplot(2,3,1);
area(dates, w_groups(squeeze(PortWS1(1,:,:)), [CORE' NEU' PIIS' EAST']));
title('GMV with sample moments')
ylim([0 1]);
subplot(2,3,4);
area(dates, w_groups(squeeze(PortWS1(2,:,:)), [CORE' NEU' PIIS' EAST']));
title('MS with sample moments')
ylim([0 1]);
subplot(2,3,2);
area(dates, w_groups(squeeze(PortWE1(1,:,:)), [CORE' NEU' PIIS' EAST']));
title('GMV with EWMA')
ylim([0 1]);
subplot(2,3,5);
area(dates, w_groups(squeeze(PortWE1(2,:,:)), [CORE' NEU' PIIS' EAST']));
title('MS with EWMA')
ylim([0 1]);
subplot(2,3,3);
area(dates, w_groups(squeeze(PortWQ1(1,:,:)), [CORE' NEU' PIIS' EAST']));
title('GMV with Equilibrium')
ylim([0 1]);
subplot(2,3,6);
area(dates, w_groups(squeeze(PortWQ1(2,:,:)), [CORE' NEU' PIIS' EAST']));
title('MS with Equilibrium')
ylim([0 1]);
legend({'CORE', 'NEU', 'PIIS', 'EAST'})
sgtitle('Group weights | No-shorting');

%% compare returns of optimal portfolio with the returns of the benchmark
CRreturns=cumprod(returns(w+1:end,:)/100+1);
% cumulated returns of the market
CRrBench=cumprod(benchReturns(w+1:end,1)/100+1);
CRS1GMV=cumprod(PortRetS1(:,1)/100+1);
CRS1MS=cumprod(PortRetS1(:,2)/100+1);
CRE1GMV=cumprod(PortRetE1(:,1)/100+1);
CRE1MS=cumprod(PortRetE1(:,2)/100+1);
CRQ1GMV=cumprod(PortRetQ1(:,1)/100+1);
CRQ1MS=cumprod(PortRetQ1(:,2)/100+1);

% GMV and MS baseline
figure;
plot(dates, [CRS1GMV CRS1MS CRrBench CREW]); %cumulated returns of sample moments compared with the market
set(gca, 'XTick', fineDateRange);
h = gca; h.XAxis.TickLabelFormat = 'M/yyyy';
grid on;
legend({'GMV','MS','Benchmark','EW'}, 'location', 'southeast');
title('Sample moments')

figure;
plot(dates, [CRE1GMV CRE1MS CRrBench CREW]); %cumulated returns of sample moments compared with the market
set(gca, 'XTick', fineDateRange);
h = gca; h.XAxis.TickLabelFormat = 'M/yyyy';
grid on;
legend({'GMV','MS','Benchmark','EW'}, 'location', 'southeast');
title('EWMA')

figure;
plot(dates, [CRQ1GMV CRQ1MS CRrBench CREW]); %cumulated returns of sample moments compared with the market
set(gca, 'XTick', fineDateRange);
h = gca; h.XAxis.TickLabelFormat = 'M/yyyy';
grid on;
legend({'GMV','MS','Benchmark','EW'}, 'location', 'southeast');
title('Equilibrium')

%%
figure
subplot(3,1,1);
plot(dates, [CRS1GMV CRS1MS CRrBench CREW]); %cumulated returns of sample moments compared with the market
%set(gca, 'XTick', fineDateRange);
%h = gca; h.XAxis.TickLabelFormat = 'M/yyyy';
grid on;
title('Sample moments')
subplot(3,1,2);
plot(dates, [CRE1GMV CRE1MS CRrBench CREW]); %cumulated returns of sample moments compared with the market
%set(gca, 'XTick', fineDateRange);
%h = gca; h.XAxis.TickLabelFormat = 'M/yyyy';
grid on;
title('EWMA')
subplot(3,1,3);
plot(dates, [CRQ1GMV CRQ1MS CRrBench CREW]); %cumulated returns of sample moments compared with the market
%set(gca, 'XTick', fineDateRange);
%h = gca; h.XAxis.TickLabelFormat = 'M/yyyy';
grid on;
title('Equilibrium')
sgtitle("Cumulated returns [%] | No-shorting");
legend({'GMV','MS','Bench','EW'}, 'location', 'eastoutside');

%% CASE 2 UPPER and/or LOWER BOUND

lowerBound = 0;
upperBound = 0.1;
lbubString = [' lb ' num2str(lowerBound*100) '% ub ' num2str(upperBound*100) '%'];

% SAMPLE MOMENTS
pp=Portfolio(p,'LowerBound',lowerBound);
pp=Portfolio(pp,'UpperBound',upperBound);
PortRetS2=zeros(r-w,2); %returns
PortWS2=zeros(2,r-w,c); %weights
pwgtMS=zeros(c,1);
% EWMA
p2=Portfolio(p1,'LowerBound',lowerBound);
p2=Portfolio(p2,'UpperBound',upperBound);
PortRetE2=zeros(r-w,2); %returns
PortWE2=zeros(2,r-w,c); %weights
pwgtMS1=zeros(c,1);
% Equilibrium
p3=Portfolio(p1,'LowerBound',lowerBound);
p3=Portfolio(p3,'UpperBound',upperBound);
PortRetQ2=zeros(r-w,2); %returns
PortWQ2=zeros(2,r-w,c); %weights
pwgtMS2=zeros(c,1);
for j=1:r-w
    % SAMPLE MOMENTS
    pp=pp.setAssetMoments(ErS(j,:),squeeze(EvS(j,:,:)));
    % GMV
    [PortWS2(1,j,:), PortRetS2(j,1)]=funGMV(pp,returns,j,w);
    % MaxSharpe
    [PortWS2(2,j,:), PortRetS2(j,2), pwgtMS]=funMS(pp,returns,j,w,pwgtMS,ErS);
    
    % EQUILIBRIUM RETURNS AND EWMA COVARIANCES
    p2=p2.setAssetMoments(ErE(j,:),squeeze(EvW(j,:,:)));
    % GMV
    [PortWE2(1,j,:), PortRetE2(j,1)]=funGMV(p2,returns,j,w);
    % MaxSharpe
    [PortWE2(2,j,:), PortRetE2(j,2), pwgtMS1]=funMS(p2,returns,j,w,pwgtMS1,ErE);
    
    % EQUILIBRIUM RETURNS AND COVARIANCES
    p3=p3.setAssetMoments(ErE(j,:),squeeze(EvE(j,:,:)));
    % GMV
    [PortWQ2(1,j,:), PortRetQ2(j,1)]=funGMV(p3,returns,j,w);
    % MaxSharpe
    [PortWQ2(2,j,:), PortRetQ2(j,2), pwgtMS2]=funMS(p3,returns,j,w,pwgtMS2,ErE);
end

%% area plots for weights over time case 2
figure;
subplot(2,3,1);
area(dates, squeeze(PortWS2(1,:,:)));
title('GMV with sample moments')
ylim([0 1]);
subplot(2,3,4);
area(dates, squeeze(PortWS2(2,:,:)));
title('MS with sample moments')
ylim([0 1]);
subplot(2,3,2);
area(dates, squeeze(PortWE2(1,:,:)));
title('GMV with EWMA')
ylim([0 1]);
subplot(2,3,5);
area(dates, squeeze(PortWE2(2,:,:)));
title('MS with EWMA')
ylim([0 1]);
subplot(2,3,3);
area(dates, squeeze(PortWQ2(1,:,:)));
title('GMV with Equilibrium')
ylim([0 1]);
subplot(2,3,6);
area(dates, squeeze(PortWQ2(2,:,:)));
title('MS with Equilibrium')
ylim([0 1]);
legend(labels)
sgtitle(['Weights |' lbubString]);

% aggregated group weights
figure;
subplot(2,3,1);
area(dates, w_groups(squeeze(PortWS2(1,:,:)), [CORE' NEU' PIIS' EAST']));
title('GMV with sample moments')
ylim([0 1]);
subplot(2,3,4);
area(dates, w_groups(squeeze(PortWS2(2,:,:)), [CORE' NEU' PIIS' EAST']));
title('MS with sample moments')
ylim([0 1]);
subplot(2,3,2);
area(dates, w_groups(squeeze(PortWE2(1,:,:)), [CORE' NEU' PIIS' EAST']));
title('GMV with EWMA')
ylim([0 1]);
subplot(2,3,5);
area(dates, w_groups(squeeze(PortWE2(2,:,:)), [CORE' NEU' PIIS' EAST']));
title('MS with EWMA')
ylim([0 1]);
subplot(2,3,3);
area(dates, w_groups(squeeze(PortWQ2(1,:,:)), [CORE' NEU' PIIS' EAST']));
title('GMV with Equilibrium')
ylim([0 1]);
subplot(2,3,6);
area(dates, w_groups(squeeze(PortWQ2(2,:,:)), [CORE' NEU' PIIS' EAST']));
title('MS with Equilibrium')
ylim([0 1]);
legend({'CORE', 'NEU', 'PIIS', 'EAST'})
sgtitle(['Group weights |' lbubString]);


%% compare returns of optimal portfolio with the returns of the benchmark
CRS2GMV=cumprod(PortRetS2(:,1)/100+1);
CRS2MS=cumprod(PortRetS2(:,2)/100+1);
CRE2GMV=cumprod(PortRetE2(:,1)/100+1);
CRE2MS=cumprod(PortRetE2(:,2)/100+1);
CRQ2GMV=cumprod(PortRetQ2(:,1)/100+1);
CRQ2MS=cumprod(PortRetQ2(:,2)/100+1);

% plots
% GMV and MS baseline
figure;
plot(dates, [CRS2GMV CRS2MS CRrBench CREW]); %cumulated returns of sample moments compared with the market
set(gca, 'XTick', fineDateRange);
h = gca; h.XAxis.TickLabelFormat = 'M/yyyy';
grid on;
legend({'GMV','MS','Benchmark','EW'}, 'location', 'southeast');
title(['Sample moments' lbubString])

figure;
plot(dates, [CRE2GMV CRE2MS CRrBench CREW]); %cumulated returns of sample moments compared with the market
set(gca, 'XTick', fineDateRange);
h = gca; h.XAxis.TickLabelFormat = 'M/yyyy';
grid on;
legend({'GMV','MS','Benchmark','EW'}, 'location', 'southeast');
title(['EWMA' lbubString])

%%
figure
subplot(3,1,1);
plot(dates, [CRS2GMV CRS2MS CRrBench CREW]); %cumulated returns of sample moments compared with the market
set(gca, 'XTick', fineDateRange);
h = gca; h.XAxis.TickLabelFormat = 'M/yyyy';
grid on;
title('Sample moments')
subplot(3,1,2);
plot(dates, [CRE2GMV CRE2MS CRrBench CREW]); %cumulated returns of sample moments compared with the market
set(gca, 'XTick', fineDateRange);
h = gca; h.XAxis.TickLabelFormat = 'M/yyyy';
grid on;
title('EWMA')
subplot(3,1,3);
plot(dates, [CRQ2GMV CRQ2MS CRrBench CREW]); %cumulated returns of sample moments compared with the market
set(gca, 'XTick', fineDateRange);
h = gca; h.XAxis.TickLabelFormat = 'M/yyyy';
grid on;
title('Equilibrium')
sgtitle(['Cumulated returns [%] | ' lbubString]);
legend({'GMV','MS','Bench','EW'}, 'location', 'eastoutside');

%% CASE 3 GROUP CONSTRAINTS

groupsBounds = [EAST; CORE; PIIS; NEU];
groupsLBound = [0.2 0.5 0 0.1];
groupsUBound = [0.5 1 0.2 0.5];

% SAMPLE MOMENTS
p=setGroups(p1,groupsBounds,groupsLBound,groupsUBound);
PortRetS3=zeros(r-w,2); %returns
PortWS3=zeros(2,r-w,c); %weights
pwgtMS=zeros(c,1);
% EWMA
p2=setGroups(p1,groupsBounds,groupsLBound,groupsUBound);
PortRetE3=zeros(r-w,2); %returns
PortWE3=zeros(2,r-w,c); %weights
pwgtMS1=zeros(c,1);
% Equilibrium
p3=setGroups(p1,groupsBounds,groupsLBound,groupsUBound);
PortRetQ3=zeros(r-w,2); %returns
PortWQ3=zeros(2,r-w,c); %weights
pwgtMS2=zeros(c,1);
for j=1:r-w
    % SAMPLE MOMENTS
    p=p.setAssetMoments(ErS(j,:),squeeze(EvS(j,:,:)));
    % GMV
    [PortWS3(1,j,:), PortRetS3(j,1)]=funGMV(p,returns,j,w);
    % MaxSharpe
    [PortWS3(2,j,:), PortRetS3(j,2), pwgtMS]=funMS(p,returns,j,w,pwgtMS,ErS);
    
    % EQUILIBRIUM RETURNS AND EWMA COVARIANCES
    p2=p2.setAssetMoments(ErE(j,:),squeeze(EvW(j,:,:)));
    % GMV
    [PortWE3(1,j,:), PortRetE3(j,1)]=funGMV(p2,returns,j,w);
    % MaxSharpe
    [PortWE3(2,j,:), PortRetE3(j,2), pwgtMS1]=funMS(p2,returns,j,w,pwgtMS1,ErE);
    
    % EQUILIBRIUM RETURNS AND COVARIANCES
    p3=p3.setAssetMoments(ErE(j,:),squeeze(EvE(j,:,:)));
    % GMV
    [PortWQ3(1,j,:), PortRetQ3(j,1)]=funGMV(p3,returns,j,w);
    % MaxSharpe
    [PortWQ3(2,j,:), PortRetQ3(j,2), pwgtMS2]=funMS(p3,returns,j,w,pwgtMS2,ErE);
end

%% area plots for weights over time case 2
% assets weights
figure;
subplot(2,3,1);
area(dates, squeeze(PortWS3(1,:,:)));
title('GMV with sample moments')
ylim([0 1]);
subplot(2,3,4);
area(dates, squeeze(PortWS3(2,:,:)));
title('MS with sample moments')
ylim([0 1]);
subplot(2,3,2);
area(dates, squeeze(PortWE3(1,:,:)));
title('GMV with EWMA')
ylim([0 1]);
subplot(2,3,5);
area(dates, squeeze(PortWE3(2,:,:)));
title('MS with EWMA')
ylim([0 1]);
subplot(2,3,3);
area(dates, squeeze(PortWQ3(1,:,:)));
title('GMV with Equilibrium')
ylim([0 1]);
subplot(2,3,6);
area(dates, squeeze(PortWQ3(2,:,:)));
title('MS with Equilibrium')
ylim([0 1]);
legend(labels)
sgtitle('Assets weights | EAST[0.2-0.5], CORE[0.5-1], PIIS[0-0.2], NEU[0.1-0.5]');

% aggregated group weights
figure;
subplot(2,3,1);
area(dates, w_groups(squeeze(PortWS3(1,:,:)), [CORE' NEU' PIIS' EAST']));
title('GMV with sample moments')
ylim([0 1]);
subplot(2,3,4);
area(dates, w_groups(squeeze(PortWS3(2,:,:)), [CORE' NEU' PIIS' EAST']));
title('MS with sample moments')
ylim([0 1]);
subplot(2,3,2);
area(dates, w_groups(squeeze(PortWE3(1,:,:)), [CORE' NEU' PIIS' EAST']));
title('GMV with EWMA')
ylim([0 1]);
subplot(2,3,5);
area(dates, w_groups(squeeze(PortWE3(2,:,:)), [CORE' NEU' PIIS' EAST']));
title('MS with EWMA')
ylim([0 1]);
subplot(2,3,3);
area(dates, w_groups(squeeze(PortWQ3(1,:,:)), [CORE' NEU' PIIS' EAST']));
title('GMV with Equilibrium')
ylim([0 1]);
subplot(2,3,6);
area(dates, w_groups(squeeze(PortWQ3(2,:,:)), [CORE' NEU' PIIS' EAST']));
title('MS with Equilibrium')
ylim([0 1]);
legend({'CORE', 'NEU', 'PIIS', 'EAST'})
sgtitle('Group weights | EAST[0.2-0.5], CORE[0.5-1], PIIS[0-0.2], NEU[0.1-0.5]');


%% compare returns of optimal portfolio with the returns of the benchmark
CRS3GMV=cumprod(PortRetS3(:,1)/100+1);
CRS3MS=cumprod(PortRetS3(:,2)/100+1);
CRE3GMV=cumprod(PortRetE3(:,1)/100+1);
CRE3MS=cumprod(PortRetE3(:,2)/100+1);
CRQ3GMV=cumprod(PortRetQ3(:,1)/100+1);
CRQ3MS=cumprod(PortRetQ3(:,2)/100+1);

%%
figure
subplot(3,1,1);
plot(dates, [CRS3GMV CRS3MS CRrBench CREW]); %cumulated returns of sample moments compared with the market
set(gca, 'XTick', fineDateRange);
h = gca; h.XAxis.TickLabelFormat = 'M/yyyy';
grid on;
title('Sample moments')
subplot(3,1,2);
plot(dates, [CRE3GMV CRE3MS CRrBench CREW]); %cumulated returns of sample moments compared with the market
set(gca, 'XTick', fineDateRange);
h = gca; h.XAxis.TickLabelFormat = 'M/yyyy';
grid on;
title('EWMA')
subplot(3,1,3);
plot(dates, [CRQ3GMV CRQ3MS CRrBench CREW]); %cumulated returns of sample moments compared with the market
set(gca, 'XTick', fineDateRange);
h = gca; h.XAxis.TickLabelFormat = 'M/yyyy';
grid on;
title('Equilibrium')
sgtitle('Cumulated returns [%] | EAST[0.2-0.5], CORE[0.5-1], PIIS[0-0.2], NEU[0.1-0.5]');
legend({'GMV','MS','Bench','EW'}, 'location', 'eastoutside');


%% Save data for further analysis
save analysis_a003.mat...
    PortRetEW...
    PortWS1 PortWE1 PortWQ1 PortWS2 PortWE2 PortWQ2 PortWS3 PortWE3 PortWQ3...
    PortRetS1 PortRetE1 PortRetQ1 PortRetS2 PortRetE2 PortRetQ2 PortRetS3 PortRetE3 PortRetQ3...
    analysisWindow w;
