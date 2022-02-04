%% 5) Advanced asset allocation
clear
clc
load data.mat

w=60;   % 5 years

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

% Daily returns (used for realized risk contribution)
returnsDaily = computeReturns(dm.dailyPrices);
datesDaily = dm.dailyDates(2:end);

% EW
PortRetEW=sum(returns(w+1:r,:),2)./c;
CREW=cumprod(PortRetEW(:,1)/100+1);
% benchmark
CRrBench=cumprod(benchReturns(w+1:end,1)/100+1);

%% Variabile initialization

% init sample moments
ErS=zeros(r-w,c);
EvS=zeros(r-w,c,c);
% constant correlation model with/without groups
EvCC2=zeros(r-w,c,c);
EvCC=zeros(r-w,c,c);

% init ex-post RC (realized)
EvSdaily=zeros(r-w,c,c);
% absolute
PortRCpost=zeros(2,r-w,c); % 2 for ERC and GRC
% relative
PortRCrpost=zeros(2,r-w,c);
% init ex-ante RC (expected)
PortRC=zeros(2,r-w,c);
PortRCr=zeros(2,r-w,c);

% init port returns and weights
% Risk Contribution
PortRetRC=zeros(r-w,2);
PortWRC=zeros(2,r-w,c);
% Criterion Function
PortRetCF=zeros(r-w,1);
PortWCF=zeros(1,r-w,c);

% starting weights
x0=ones(c,3)/c; % 3 because we have 3 strategies using fmincon
ew=ones(c,1)/c;

% optimizazion options
options = optimoptions(@fmincon,'MaxIterations',3000,...
                       'MaxFunctionEvaluations',10000,'Display','off');
optionsGRC=options;

%% Constraints
% fully invested, w sum to 1
Aeq=ones(1,c);
beq=1;

% lower/upper bounds
lb=ones(c,1)*(0.0);
ub=ones(c,1)*(1);

% Four groups for RC:
% 1) CZ, HN, PO      2) ES, IT, IR, PT
% 3) UK, NW, SW, SD  4) the rest
global Gmat;
PIIS = zeros(18,1)'; PIIS(labels == "PT" | labels == "IT" | labels == "ES" | labels == "IR") = 1;
EAST = zeros(18,1)'; EAST(labels == "CZ" | labels == "PO" | labels == "HN") = 1;
NEU = zeros(18,1)'; NEU(labels == "UK" | labels == "NW" | labels == "SD" | labels == "SW") = 1;
CORE = ones(18,1)' - EAST - PIIS - NEU;
Gmat = [EAST; CORE; PIIS; NEU];
% group risk budgets
global bvec;
% east countries 20% of total risk, core 50% of tot risk...
bvec = [0.2 0.5 0.1 0.2];

%% Inputs
% rolling sample moments
for j=1:r-w
    curReturns = returns(j:j+w-1,:);
    % Monthly measures (for EX-ANTE computation)
    % sample mean
    ErS(j,:)=mean(curReturns);
    % sample varcov
    EvS(j,:,:)=cov(curReturns);
    % CCM with groups
    EvCC2(j,:,:) = const_corr2(returns(j:j+w-1,:),[CORE' PIIS' NEU' EAST']);
    % CCM 
    EvCC(j,:,:) = const_corr2(returns(j:j+w-1,:));
    
    % Daily measures (for EX-POST computation)    
    % ErS(j,) and EvS(j,) are the estimates for period date(j).
    % Calculate last day of the month we are doing the estimation for
    finalDate = eomdate(dates(j));
    % calculate first date as 1 months before finalDate
    initialDate = finalDate - calmonths(1) + caldays(1);
    % restrict daily prices range and calculate varcov matrix
    dataFilter = (datesDaily >= initialDate) & (datesDaily <= finalDate);
    returnsRg = returnsDaily(dataFilter, :);
    EvSdaily(j,:,:) = cov(returnsRg);
    % recall that EvSdaily are not used as estimates for next period
end

% check that matrices are positive definite
cSample = 0;
cDaily = 0;
cCC = 0;
for j=1:size(EvS,1)
    cSample = all(eig(squeeze(EvS(j,:,:)))>eps) + cSample;
    cDaily = all(eig(squeeze(EvSdaily(j,:,:)))>eps) + cDaily;
    cCC = all(eig(squeeze(EvCC(j,:,:)))>eps) + cCC;
end
if (cSample ~= 60 || cDaily ~= 60 || cCC ~= 60)
    disp('Some varcov matrix is not positive definite')
end
% all matrices are positive definite
clear cSample cDaily cCC


% RNN estimated returns
nnReturns = csvread('returnPredsRNN.csv');
PortWNN = csvread('weightsRNN.csv');
lastRet = returns(end-w+1:end,:);
PortRetNN = sum(lastRet.*PortWNN,2);
CRNN=cumprod(PortRetNN/100+1);
clear lastRet;

%% Optimization routine
optimTime = 0;
for j=1:r-w
    disp(['Iteration:' num2str(j)])
    % batches of 60 returns up to r-1, used to estimate weights of next
    % period (up to r)
    curReturns = returns(j:j+w-1,:);
    curBenchReturns = benchReturns(j:j+w-1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [1] custom criterion function
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % strategy identifier
    currStr=1;
    % maxdd and ES uses median instead of mean to reduce turnover
    critW=@(x)(-expshortfall(x,curReturns,0.05)-maxdd(x,curReturns));
    
    tic;
    wCrit=fmincon(critW,x0(:,currStr),[],[],Aeq,beq,lb,ub,[],options);
    optimTime=optimTime+toc;
    % store weights and returns
    PortWCF(1,j,:)=wCrit';
    PortRetCF(j,1)=returns(j+w,:)*wCrit;
    % store current weights, used as x0 of next iteration
    x0(:,currStr) = wCrit;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [2] equal risk contribution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    currStr=2;
    global Sigma;
    Sigma=squeeze(EvS(j,:,:));
    ercW = @(w) ERCfun(w,Sigma); 
    tic;
    wERC = fmincon(ercW,x0(:,currStr),[],[],Aeq,beq,lb,ub,[],options);
    optimTime=optimTime+toc;
    % store weights and returns
    PortWRC(1,j,:)=wERC'/sum(wERC);
    PortRetRC(j,1)=returns(j+w,:)*wERC;
    % store current weights, used as x0 of next iteration
    x0(:,currStr) = wERC;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [3] general risk contribution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    currStr=3;
    tic;
    wGRC = fmincon(critW,x0(:,currStr),[],[],Aeq,beq,lb,ub,@GRCcon,optionsGRC);
    optimTime=optimTime+toc;
    % store weights and returns
    PortWRC(2,j,:)=wGRC'/sum(wGRC);
    PortRetRC(j,2)=returns(j+w,:)*wGRC;
    % store current weights, used as x0 of next iteration
    x0(:,currStr) = wGRC;
    
    % check estimated optimization time
    cycles = r-w-1;
    totalDuration = round(sum(optimTime)*cycles);
    if j == 1
        question = ['The optimization process could take up to ', num2str(totalDuration), 's. Continue?'];
        choice = menu(question,'Yes','No');
        if choice==2
            disp('Optimization process aborted by user');
            break;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Risk contributions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EX-ANTE (expected risk contributions)
    % ERC
    RiskMeasure = squeeze(EvS(j,:,:));
    PortRC(1,j,:)=(1/sqrt(wERC'*RiskMeasure*wERC))*((RiskMeasure*wERC).*wERC);
    PortRCr(1,j,:)=PortRC(1,j,:)/sum(PortRC(1,j,:));
    % GRC
    PortRC(2,j,:)=(1/sqrt(wGRC'*RiskMeasure*wGRC))*((RiskMeasure*wGRC).*wGRC);
    PortRCr(2,j,:)=PortRC(2,j,:)/sum(PortRC(2,j,:));

    % EX-POST (realized risk contributions)
    RiskMeasure = squeeze(EvSdaily(j,:,:));
    % ERC
    PortRCpost(1,j,:)=(1/sqrt(wERC'*RiskMeasure*wERC))*((RiskMeasure*wERC).*wERC);
    PortRCrpost(1,j,:)=PortRCpost(1,j,:)/sum(PortRCpost(1,j,:));
    % GRC
    PortRCpost(2,j,:)=(1/sqrt(wGRC'*RiskMeasure*wGRC))*((RiskMeasure*wGRC).*wGRC);
    PortRCrpost(2,j,:)=PortRCpost(2,j,:)/sum(PortRCpost(2,j,:));
    
    if any(corrcov(RiskMeasure) < 0, 'all')
        disp(['Negative correlation found at iteration ' num2str(j)])
    end
    
    if any(wERC<0)
        disp('At least one weight of wERC is negative');
        break;
    end
end

%% Other
p=Portfolio;
p=p.setAssetList(labels);
p=p.setDefaultConstraints;

% NN
PortRetNN1=zeros(r-w,2);
PortWNN1=zeros(2,r-w,c);
pwgtMSNN1=zeros(c,1);
lowerBound = 0;
upperBound = 0.1;
p1=Portfolio(p,'LowerBound',lowerBound);
p1=Portfolio(p1,'UpperBound',upperBound);

% CCM 
PortRetCC=zeros(r-w,2);
PortWCC=zeros(2,r-w,c);
pwgtMSCC=zeros(c,1);
lowerBound = 0;
upperBound = 0.1;
p2=Portfolio(p,'LowerBound',lowerBound);
p2=Portfolio(p2,'UpperBound',upperBound);

% CCM with groups 
PortRetCC2=zeros(r-w,2);
PortWCC2=zeros(2,r-w,c);
pwgtMSCC2=zeros(c,1);
lowerBound = 0;
upperBound = 0.1;
p3=Portfolio(p,'LowerBound',lowerBound);
p3=Portfolio(p3,'UpperBound',upperBound);

for j=1:r-w
    % NN
    p1=p1.setAssetMoments(nnReturns(j,:),squeeze(EvS(j,:,:)));
    % GMV
    [PortWNN1(1,j,:), PortRetNN1(j,1)]=funGMV(p1,returns,j,w);
    % MaxSharpe
    [PortWNN1(2,j,:), PortRetNN1(j,2), pwgtMSNN1]=funMS(p1,returns,j,w,pwgtMSNN1,nnReturns);
    
    % CCM with groups
    p3=p3.setAssetMoments(ErS(j,:),squeeze(EvCC2(j,:,:)));
    % GMV
    [PortWCC2(1,j,:), PortRetCC2(j,1)]=funGMV(p3,returns,j,w);
    % MS
    [PortWCC2(2,j,:), PortRetCC2(j,2),pwgtMSCC]=funMS(p3,returns,j,w,pwgtMSCC,ErS);
    
    % CCM 
    p2=p2.setAssetMoments(ErS(j,:),squeeze(EvCC(j,:,:)));
    % GMV
    [PortWCC(1,j,:), PortRetCC(j,1)]=funGMV(p2,returns,j,w);
    % MS
    [PortWCC(2,j,:), PortRetCC(j,2),pwgtMSCC]=funMS(p2,returns,j,w,pwgtMSCC,ErS);
end

%% RC over time
figure
subplot(2,2,1)
area(dates,squeeze(PortRCr(1,:,:)))
datetick('x', 'yyyy')
title('EX-ANTE ERC')
subplot(2,2,2)
W=squeeze(PortRCrpost(1,:,:));
area(dates,W.*(W>0));
hold on;
area(dates,W.*(W<0));
datetick('x', 'yyyy')
title('EX-POST ERC')
hold off;
subplot(2,2,3)
area(dates,squeeze(PortRCr(2,:,:)))
datetick('x', 'yyyy')
title('EX-ANTE GRC (EAST:20%, CORE:50%, PIIS:10%, NEU:20%)')
subplot(2,2,4)
W=squeeze(PortRCrpost(2,:,:));
area(dates,W.*(W>0));
hold on;
area(dates,W.*(W<0));
datetick('x', 'yyyy')
title('EX-POST GRC (EAST:20%, CORE:50%, PIIS:10%, NEU:20%)')
legend(labels)
hold off;
sgtitle("Risk contribution over time");

%% Weights over time
figure
subplot(2,2,1)
area(dates,squeeze(PortWCF(1,:,:)))
title('Custom criterion')
ylim([0 1]);
subplot(2,2,2)
area(dates,squeeze(PortWRC(1,:,:)))
title('ERC')
ylim([0 1]);
subplot(2,2,3)
area(dates,squeeze(PortWRC(2,:,:)))
title('GRC (EAST:20%, CORE:50%, PIIS:10%, NEU:20%)')
ylim([0 1]);
subplot(2,2,4)
area(dates,PortWNN)
legend(labels)
title('RNN')
ylim([0 1]);
sgtitle("Weights over time");

% aggregated group weights
figure;
subplot(2,2,1);
area(dates, w_groups(squeeze(PortWCF(1,:,:)), [CORE' NEU' PIIS' EAST']));
title('Custom criterion')
ylim([0 1]);
subplot(2,2,2);
area(dates, w_groups(squeeze(PortWRC(1,:,:)), [CORE' NEU' PIIS' EAST']));
title('ERC')
ylim([0 1]);
subplot(2,2,3);
area(dates, w_groups(squeeze(PortWRC(2,:,:)), [CORE' NEU' PIIS' EAST']));
title('GRC (EAST:20%, CORE:50%, PIIS:10%, NEU:20%)')
ylim([0 1]);
subplot(2,2,4);
area(dates, w_groups(PortWNN, [CORE' NEU' PIIS' EAST']));
title('RNN')
ylim([0 1]);
legend({'CORE', 'NEU', 'PIIS', 'EAST'})
sgtitle('Group weights over time');

% weights
figure
subplot(2,3,1)
area(dates,squeeze(PortWNN1(1,:,:)))
title('RNN GMV')
ylim([0 1]);
subplot(2,3,4)
area(dates,squeeze(PortWNN1(2,:,:)))
title('RNN MS')
ylim([0 1]);
subplot(2,3,2)
area(dates,squeeze(PortWCC2(1,:,:)))
title('CCM with groups GMV')
ylim([0 1]);
subplot(2,3,5)
area(dates,squeeze(PortWCC2(2,:,:)))
title('CCM with groups MS')
ylim([0 1]);
subplot(2,3,3)
area(dates,squeeze(PortWCC(1,:,:)))
title('CCM GMV')
ylim([0 1]);
subplot(2,3,6)
area(dates,squeeze(PortWCC(2,:,:)))
title('CCM MS')
ylim([0 1]);
legend(labels)
sgtitle("Weights over time");

% aggregated group weights
figure;
subplot(2,3,1);
area(dates, w_groups(squeeze(PortWNN1(1,:,:)), [CORE' NEU' PIIS' EAST']));
title('RNN GMV')
ylim([0 1]);
subplot(2,3,4);
area(dates, w_groups(squeeze(PortWNN1(2,:,:)), [CORE' NEU' PIIS' EAST']));
title('RNN MS')
ylim([0 1]);
subplot(2,3,2);
area(dates, w_groups(squeeze(PortWCC2(1,:,:)), [CORE' NEU' PIIS' EAST']));
title('CCM with groups GMV')
ylim([0 1]);
subplot(2,3,5);
area(dates, w_groups(squeeze(PortWCC2(2,:,:)), [CORE' NEU' PIIS' EAST']));
title('CCM with groups MS')
ylim([0 1]);
subplot(2,3,3);
area(dates, w_groups(squeeze(PortWCC(1,:,:)), [CORE' NEU' PIIS' EAST']));
title('CCM GMV')
ylim([0 1]);
subplot(2,3,6);
area(dates, w_groups(squeeze(PortWCC(2,:,:)), [CORE' NEU' PIIS' EAST']));
title('CCM MS')
ylim([0 1]);
legend({'CORE', 'NEU', 'PIIS', 'EAST'})
sgtitle('Group weights over time');

%% Port cumulated returns
% gmv and ms of CCM, CCM2 and NN1
CRNN1GMV=cumprod(PortRetNN1(:,1)/100+1);
CRNN1MS=cumprod(PortRetNN1(:,2)/100+1);
CRCCGMV=cumprod(PortRetCC(:,1)/100+1);
CRCCMS=cumprod(PortRetCC(:,2)/100+1);
CRCC2GMV=cumprod(PortRetCC2(:,1)/100+1);
CRCC2MS=cumprod(PortRetCC2(:,2)/100+1);

figure
subplot(3,1,1);
plot(dates, [CRNN1GMV CRNN1MS CRrBench CREW]);
set(gca, 'XTick', fineDateRange);
h = gca; h.XAxis.TickLabelFormat = 'M/yyyy';
grid on;
title('NN1')
subplot(3,1,2);
plot(dates, [CRCCGMV CRCCMS CRrBench CREW]);
set(gca, 'XTick', fineDateRange);
h = gca; h.XAxis.TickLabelFormat = 'M/yyyy';
grid on;
title('CCM')
subplot(3,1,3);
plot(dates, [CRCC2GMV CRCC2MS CRrBench CREW]);
set(gca, 'XTick', fineDateRange);
h = gca; h.XAxis.TickLabelFormat = 'M/yyyy';
grid on;
title('CCM with groups')
sgtitle("Cumulated returns [%]");
legend({'GMV','MS','Bench','EW'}, 'location', 'eastoutside');

% CR criterion function, RC, NN
CRCF=cumprod(PortRetCF(:,1)/100+1);
CRERC=cumprod(PortRetRC(:,1)/100+1);
CRGRC=cumprod(PortRetRC(:,2)/100+1);

figure
plot(dates, [CRCF CRERC CRGRC CRNN CRrBench CREW]);
set(gca, 'XTick', fineDateRange);
h = gca; h.XAxis.TickLabelFormat = 'M/yyyy';
grid on;
title('Cumulated returns [%]')
legend({'CF','ERC','GRC','NN','Bench','EW'}, 'location', 'eastoutside');

%% Save data for further analysis
save analysis_a004.mat...
    PortRetEW...
    PortWCF PortWRC PortWNN PortWNN1 PortWCC2 PortWCC...
    PortRetCF PortRetRC PortRetNN PortRetNN1 PortRetCC2 PortRetCC...
    analysisWindow w;