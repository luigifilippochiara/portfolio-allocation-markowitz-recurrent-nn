%% Track record analysis
clear
clc
load data.mat
load analysis_a004.mat

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
% CF = Custom criterion function
% ERC = Equal risk contribution
% GRC = General risk contribution (minimizing CF, with group risk budgets)
% CC = Constant correlation model
% CC2 = Constant correlation with groups model
% NN1 = Recurrent neural network in Markovitz framework
% NN = Recurrent neural network

% Strategy names
strategyNames = {'CF', 'ERC', 'GRC',...
                 'CC', 'CC2', 'NN1',...
                 'NN'};

benchReturns5yrs = benchReturns(end-analysisWindow+1:end);

allPortW = {...
    squeeze(PortWCF(1,:,:)) squeeze(PortWRC(1,:,:)) squeeze(PortWRC(2,:,:)) squeeze(PortWCC(1,:,:)) squeeze(PortWCC(2,:,:))...
    squeeze(PortWCC2(1,:,:)) squeeze(PortWCC2(2,:,:)) squeeze(PortWNN1(1,:,:)) squeeze(PortWNN1(2,:,:)),...
    PortWNN};

%% Turnover
% focusing on the fraction of portfolio which is changing
tRows = size(squeeze(PortWCF(1,:,:)), 1); 
tCols = length(allPortW);
turnovers = zeros(tRows-1, tCols); % lose one row because we compute lagged difference
effTurnovers = zeros(tRows-1, tCols); % lose one row because we compute lagged difference
for i=1:length(allPortW)
    turnovers(:,i) = approxTurnover(allPortW{i});
    effTurnovers(:,i) = effectiveTurnover(allPortW{i},returns(end-w+1:end,:));
end
effTurnoverEWmean = mean(effectiveTurnover(PortWEW,returns(end-w+1:end,:)));

% plots
legendNames = {...
    strategyNames{1};strategyNames{2};strategyNames{3};['GMV ' strategyNames{4}];['MS ' strategyNames{4}];...
    ['GMV ' strategyNames{5}];['MS ' strategyNames{5}];['GMV ' strategyNames{6}];['MS ' strategyNames{6}];...
    strategyNames{7}};
figure
plot(dates(2:end), turnovers);
grid on;
legend(legendNames);
title('Turnover on all strategies');

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
plot(dates, div);
h = gca; h.XAxis.TickLabelFormat = 'M/yyyy';
grid on;
legend(legendNames);
title('Concentration index on all strategies');

divMetrics.mean = mean(div, 1);
divMetrics.max = max(div, [], 1);
divMetrics.min = min(div, [], 1);

%% Descriptive analysis of portfolios
strategyColNames = ['Benchmark';'EW';legendNames];

allRet=[benchReturns5yrs PortRetEW...
        PortRetCF(:,1) PortRetRC(:,1) PortRetRC(:,2) PortRetCC(:,1) PortRetCC(:,2)...
        PortRetCC2(:,1) PortRetCC2(:,2) PortRetNN1(:,1) PortRetNN1(:,2)...
        PortRetNN];

% Computing moments
allRetMetrics.mean=mean(allRet)';
allRetMetrics.std=sqrt(var(allRet)');
% Other quantities
allRetMetrics.min=min(allRet)';
allRetMetrics.max=max(allRet)';

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
plot(dates, DD(:,3:end));
grid on;
legend(strategyNames);
title('Drawdown sequence on all strategies');

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
% format bank % limit to 2 digits
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
format
