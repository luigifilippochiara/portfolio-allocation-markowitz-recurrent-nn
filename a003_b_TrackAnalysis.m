%% Track record analysis
clear
clc
load data.mat
load analysis_a003.mat

% restrict returns to last 5 years
returns = returns(end-analysisWindow-w+1:end,:);
benchReturns = benchReturns(end-analysisWindow-w+1:end,:);
% Restrict dates to last 5 years of data
dates = dates(end-analysisWindow+1:end);

% Detailed date range for bigger plots
fineDateRange = (dates(1) : 180 : dates(end));

[r,c]=size(returns);

PortWEW = ones(size(returns(end-w+1:end,:)))./c;

%% Initialize
% S1 = Sample moments - no short selling
% E1 = Equilibrium returns and EWMA covariances - no short selling
% Q1 = Equilibrium returns and covariances - no short selling
% S2 = Sample moments - lower bound + upper bound
% E2 = Eq returns, EWMA cov - lower bound + upper bound
% Q2 = Eq returns and cov - lower bound + upper bound
% S3 = Sample moments - group constraints
% E3 = Eq returns, EWMA cov - group constraints
% Q3 = Eq returns and cov - group constraints

% Strategy names
strategyNames = {'SM, no short', 'EQ+EWMA, no short', 'EQ, no short',...
                 'SM, bounds',  'EQ+EWMA, bounds', 'EQ, bounds',...
                 'SM, group bounds', 'EQ+EWMA, group bounds', 'EQ, group bounds'};

benchReturns5yrs = benchReturns(end-analysisWindow+1:end);

allPortW = {...
    squeeze(PortWS1(1,:,:)) squeeze(PortWS1(2,:,:)) squeeze(PortWE1(1,:,:)) squeeze(PortWE1(2,:,:)) squeeze(PortWQ1(1,:,:)) squeeze(PortWQ1(2,:,:))...
    squeeze(PortWS2(1,:,:)) squeeze(PortWS2(2,:,:)) squeeze(PortWE2(1,:,:)) squeeze(PortWE2(2,:,:)) squeeze(PortWQ2(1,:,:)) squeeze(PortWQ2(2,:,:))...
    squeeze(PortWS3(1,:,:)) squeeze(PortWS3(2,:,:)) squeeze(PortWE3(1,:,:)) squeeze(PortWE3(2,:,:)) squeeze(PortWQ3(1,:,:)) squeeze(PortWQ3(2,:,:))};

%% Turnover
% focusing on the fraction of portfolio which is changing
tRows = size(squeeze(PortWS1(1,:,:)), 1); 
tCols = length(allPortW);
turnovers = zeros(tRows-1, tCols); % lose one row because we compute lagged difference
effTurnovers = zeros(tRows-1, tCols); % lose one row because we compute lagged difference
for i=1:length(allPortW)
    turnovers(:,i) = approxTurnover(allPortW{i});
    effTurnovers(:,i) = effectiveTurnover(allPortW{i},returns(end-w+1:end,:));
end
effTurnoverEWmean = mean(effectiveTurnover(PortWEW,returns(end-w+1:end,:)));

% plots
% GMV and MS on S1, the baseline with SM
figure;
ax1 = subplot(2,2,1);
plot(dates(2:end), [turnovers(:,1) turnovers(:,2)]);
ax1.XAxis.TickLabelFormat = 'yyyy';
grid on;
legend({'GMV', 'MS'});
title('Approximate');
ax2 = subplot(2,2,2);
plot(dates(2:end), [effTurnovers(:,1) effTurnovers(:,2)]);
ax2.XAxis.TickLabelFormat = 'yyyy';
grid on;
legend({'GMV', 'MS'});
title('Effective');

% Compare different set of constraints in GMV and MS portfolios
% GMV
ax4 = subplot(2,2,3);
plot(dates(2:end), [turnovers(:,1) turnovers(:,3) turnovers(:,5) turnovers(:,7) turnovers(:,9) turnovers(:,11) turnovers(:,13) turnovers(:,15) turnovers(:,17)]);
ax4.XAxis.TickLabelFormat = 'yyyy';
grid on;
legend(strategyNames);
title('GMV on all strategies');
% MS
ax5 = subplot(2,2,4);
plot(dates(2:end), [turnovers(:,2) turnovers(:,4) turnovers(:,6) turnovers(:,8) turnovers(:,10) turnovers(:,12) turnovers(:,14) turnovers(:,16) turnovers(:,18)]);
ax5.XAxis.TickLabelFormat = 'yyyy';
grid on;
legend(strategyNames);
title('MS on all strategies');
linkaxes([ax1, ax2, ax4, ax5], 'y');
sgtitle('Turnover on different strategies');
clear ax1 ax2 ax4 ax5;

% Compute mean, min, max to have an idea of differences of port composition
% variation over time
turnoverMetrics.mean = mean(effTurnovers, 1);
turnoverMetrics.max = max(effTurnovers, [], 1);
turnoverMetrics.min = min(effTurnovers, [], 1);

%% Concentration index
% Is a proxy for diversification, the higher, the more concentrated the port
div = zeros(tRows, tCols);
for i=1:length(allPortW)
    div(:,i) = diversification(allPortW{i}, c);
end

figure;
ax1 = subplot(1,2,1);
plot(dates, [div(:,1) div(:,3) div(:,5) div(:,7) div(:,9) div(:,11) div(:,13) div(:,15) div(:,17)]);
ax1.XAxis.TickLabelFormat = 'yyyy';
grid on;
legend(strategyNames);
title('GMV');
ax2 = subplot(1,2,2);
plot(dates, [div(:,2) div(:,4) div(:,6) div(:,8) div(:,10) div(:,12) div(:,14) div(:,16) div(:,18)]);
ax2.XAxis.TickLabelFormat = 'yyyy';
grid on;
legend(strategyNames);
title('MS');
linkaxes([ax1, ax2], 'y');
sgtitle('Concentration index on all strategies');
clear ax1 ax2;

divMetrics.mean = mean(div, 1);
divMetrics.max = max(div, [], 1);
divMetrics.min = min(div, [], 1);

%% Descriptive analysis of portfolios
strategyColNames = {'Benchmark';'EW';...
    ['GMV ' strategyNames{1}];['MS ' strategyNames{1}];['GMV ' strategyNames{2}];['MS ' strategyNames{2}];['GMV ' strategyNames{3}];['MS ' strategyNames{3}];...
    ['GMV ' strategyNames{4}];['MS ' strategyNames{4}];['GMV ' strategyNames{5}];['MS ' strategyNames{5}];['GMV ' strategyNames{6}];['MS ' strategyNames{6}];...
    ['GMV ' strategyNames{7}];['MS ' strategyNames{7}];['GMV ' strategyNames{8}];['MS ' strategyNames{8}];['GMV ' strategyNames{9}];['MS ' strategyNames{9}]};

allRet=[benchReturns5yrs PortRetEW...
            PortRetS1(:,1) PortRetS1(:,2) PortRetE1(:,1) PortRetE1(:,2) PortRetQ1(:,1) PortRetQ1(:,2)...
            PortRetS2(:,1) PortRetS2(:,2) PortRetE2(:,1) PortRetE2(:,2) PortRetQ2(:,1) PortRetQ2(:,2)...
            PortRetS3(:,1) PortRetS3(:,2) PortRetE3(:,1) PortRetE3(:,2) PortRetQ3(:,1) PortRetQ3(:,2)];
% Computing moments
allRetMetrics.mean=mean(allRet)';
allRetMetrics.std=sqrt(var(allRet)');
% Other quantities
allRetMetrics.min=min(allRet)';
allRetMetrics.max=max(allRet)';
% allRetMetrics.q=quantile(allRet,[0.01 0.05 0.5 0.95 0.99])';

%% Drawdown sequence
DD=zeros(size(allRet));
for i=1:size(allRet,2)
    DD(1,i)=min(allRet(1,i)/100,0);
    for j=2:size(allRet,1)
        DD(j,i)=min(0,(1+DD(j-1,i))*(1+allRet(j,i)/100)-1);
    end
end
% maximum drawdown
maxDD = max(abs(DD))'*100;

figure;
ax1 = subplot(1,2,1);
plot(dates, [DD(:,3) DD(:,5) DD(:,7) DD(:,9) DD(:,11) DD(:,13) DD(:,15) DD(:,17) DD(:,19)]);
ax1.XAxis.TickLabelFormat = 'M-yyyy';
grid on;
legend(strategyNames);
title('GMV');
ax2 = subplot(1,2,2);
plot(dates, [DD(:,4) DD(:,6) DD(:,8) DD(:,10) DD(:,12) DD(:,14) DD(:,16) DD(:,18) DD(:,20)]);
ax2.XAxis.TickLabelFormat = 'M-yyyy';
grid on;
legend(strategyNames);
title('MS');
linkaxes([ax1, ax2], 'y');
sgtitle('Drawdown sequence on all strategies');
clear ax1 ax2;

%% Performance and Risk indicators
% Betas
perf.betas = zeros(size(allRet,2),1);
for i=1:size(allRet,2)
    perf.betas(i) = beta(allRet(:,i),benchReturns5yrs);
end

% Sharpe
perf.sharpe = mean(allRet)'./sqrt(var(allRet))';

% Sortino
perf.sortino = zeros(size(perf.sharpe));
% compute volatility among negative elements of all assets
for i=1:size(allRet,2)
    negativeReturns = allRet(allRet(:,i)<0, i);
    perf.sortino(i)=sqrt(var(negativeReturns));
end
clear negativeReturns;
perf.sortino = mean(allRet)'./ perf.sortino;

% Treynor
perf.treynor = mean(allRet)'./ perf.betas;

% Value at Risk (currently not used)
% alpha=0.05;
% q=quantile(allRet,alpha);
% perf.VaR = abs(q)';
% clear q alpha;

% Expected Shortfall
alpha=0.05;
perf.es = zeros(size(perf.sharpe));
for i=1:size(allRet,2)
    % filter returns smaller than alpha quantile
    filter = allRet(:,i) < quantile(allRet(:,i), alpha);
    % conditional mean on filtered returns
    perf.es(i) = mean(allRet(filter, i));
end
perf.es = abs(perf.es);
clear alpha filter;

% Calmar ratio
perf.calmar = mean(allRet)'./maxDD;
% Sterling ratio
largestDD=5;
largestDDavg=zeros(size(allRet,2),1);
for j=1:size(allRet,2)
    % average of the largest DD
    [sDDj,~]=sort(abs(DD(:,j)),'descend');
    largestDDavg(j)=mean(sDDj(1:largestDD))*100;
end
perf.sterling = mean(allRet)'./largestDDavg;
clear largestDD largestDDavg sDDj;

% Risk indicators table
colNames = {'Strategy' 'meanReturn' 'StDev' 'Beta' 'maxDD' 'ES' 'avgConcIndex' 'avgTurnover'};
riskIndTab=table(strategyColNames,...
    allRetMetrics.mean,allRetMetrics.std,...
    perf.betas,maxDD,perf.es,[median(divMetrics.mean) 0 divMetrics.mean]',[effTurnoverEWmean effTurnoverEWmean turnoverMetrics.mean]',...
    'VariableNames',colNames);
format bank % limit to 2 digits
% returnsTab
f = figure;
uit = uitable(f,'Data',table2cell(riskIndTab),...
            'ColumnName',colNames,...
            'unit','normalized','Position',[0 0 1 1]);
        
%% Ranking
% Performance indicators table
[~,sharpe_r]=sort(perf.sharpe,'descend');
sharpe_r(sharpe_r)=1:length(sharpe_r);

[~,sortino_r]=sort(perf.sortino,'descend');
sortino_r(sortino_r)=1:length(sortino_r);

[~,treynor_r]=sort(perf.treynor,'descend');
treynor_r(treynor_r)=1:length(treynor_r);

[~,sterling_r]=sort(perf.sterling,'descend');
sterling_r(sterling_r)=1:length(sterling_r);

[~,turnover_r]=sort([effTurnoverEWmean effTurnoverEWmean turnoverMetrics.mean]','ascend');
turnover_r(turnover_r)=1:length(turnover_r);

[~,div_r]=sort([median(divMetrics.mean) 0 divMetrics.mean]','ascend');
div_r(div_r)=1:length(div_r);

allPMr=[sharpe_r sortino_r treynor_r sterling_r turnover_r div_r];
% sum over rows and best is the strategy with lowest nr

% compute a composite index
CIpm=allPMr*[1 1 1 1 4 4]';


varNames = {'Strategy' 'Rank' 'Sharpe' 'Sortino' 'Treynor' 'Sterling' 'avgConcIndex' 'avgTurnover'};
perfIndTab=table(strategyColNames,...
    CIpm, perf.sharpe,perf.sortino,perf.treynor,perf.sterling,...
    [median(divMetrics.mean) 0 divMetrics.mean]', [effTurnoverEWmean effTurnoverEWmean turnoverMetrics.mean]',... 
    'VariableNames',varNames);
perfIndTab = sortrows(perfIndTab, 'Rank', 'ascend');
f = figure;
uit = uitable(f,'Data',table2cell(perfIndTab),...
            'ColumnName',varNames,...
            'unit','normalized','Position',[0 0 1 1]);
        
