%% 1) Preliminary analysis
% Portfolio overview and descriptive analysis of the dataset

clear
clc
load data.mat

% store all labels together for this initial analysis
allLabels = [labels' benchLabel];
% recall that dm.dates stores dates of raw prices, while dates variable
% stores the dates of the returns (corresponding to dm.dates[2:end]
priceDates = dm.dates;
% Detailed date range for bigger plots
fineDateRange = (priceDates(1) : years(1) : priceDates(end));
fineDateRange = [fineDateRange priceDates(end)] + caldays(4);

%% Show presence of nan values in original data
% benchmark should be on first column
figure
priceGaps = ~isnan(dm.dailyPricesWithNans);
h = imagesc((priceGaps.*1)');
set(gca,'Ytick',1:length(labels)+1,'YTickLabel', [benchLabel labels']);
% xticks with unique years from jan2001
years = unique(year(dm.datesWithNans(2:end)));
yearsNr = linspace(1,size(dm.dailyPricesWithNans, 1),length(years));
set(gca,'Xtick',yearsNr,'XTickLabel',years);
xtickangle(45)
title('Nan values in original Eikon data');

%% Plot normalized prices
% compute normalized prices (divide each series by the initial price)
pricesNorm = dm.prices./dm.prices(1,:);
benchPricesNorm = dm.benchPrices./dm.benchPrices(1,:);

% reorder columns following the order of the last price row
% this will be used to plot the legend in a more understandable order
allPricesNorm = [pricesNorm benchPricesNorm];
[~, order] = sort(allPricesNorm(end,:), 'descend');
% find EU index position
benchPosition = find(order == size(allPricesNorm,2));

figure;
p=plot(priceDates, allPricesNorm(:,order));
datePositions=[0.8 0.8 0.7 0.8 0.7 0.8 0.8];
plotEvents(datePositions, 0.09);
set(gca, 'XTick', fineDateRange);
h = gca; h.XAxis.TickLabelFormat = 'yyyy';
%xlabel('Date');
ylabel('Normalized asset prices');
title('Normalized Asset Prices and Benchmark');
grid on;
% change width of benchmark line
p(benchPosition).LineWidth = 3;
l=legend(allLabels(order),'location','eastoutside');
title(l, 'Last Price[Decr]','FontWeight','normal');

%% Plot returns
allReturns = [returns benchReturns];
% reorder the assets in decreasing variance order
allVar = diag(cov(allReturns));
[~, order] = sort(allVar, 'descend');
% find EU index position
benchPosition = find(order == length(allVar));

figure;
p=plot(dates, allReturns(:,order));
datePositions=[-14 -14 -15 -14 -15 -14];
plotEvents(datePositions, 1);
set(gca, 'XTick', fineDateRange);
h = gca; h.XAxis.TickLabelFormat = 'yyyy';
%xlabel('Date');
ylabel('Returns [%]');
title('Asset monthly returns');
grid on;
p(benchPosition).LineWidth = 3;
l=legend(allLabels(order), 'location', 'eastoutside');
title(l,'Var[Decr]','FontWeight','normal');

%% Plot moments

figure;
% first moments
subplot(1,2,1);
meanReturns = mean(returns);
benchMeanReturn = mean(benchReturns);
allMeanReturns = [meanReturns benchMeanReturn];
bar(allMeanReturns);
set(gca,'Xtick',1:length(allMeanReturns),'XTickLabel',allLabels);
title('Mean monthly returns');
ylabel('Returns [%]');
grid on;

% second moments
subplot(1,2,2);
allVarCov = cov(allReturns);
volatility = sqrt(diag(allVarCov));
bar(volatility);
set(gca,'Xtick',1:length(allLabels),'XTickLabel',allLabels);
title('Asset volatilities');
ylabel('Volatility [%]');
grid on;


%% Plot asset return distributions and qqplots
% maybe could be useful to see how
% far we are from normality and if it is skewed towards positive or
% negative returns
% examples:
% ttps://medium.com/@mariano.scandizzo/strategic-asset-allocation-with-python-c9afef392e90
% https://orfe.princeton.edu/~jqfan/fan/FinEcon/chap1.pdf p.25
plotCols = 4;
plotRows = ceil(length(allMeanReturns)/plotCols);
% disp(plotRows)
for i = 1:plotRows
    for j = 1:plotCols
        plotIndex = ((i-1)*plotCols)+j;
        % terminate if last asset has been plot
        if plotIndex > length(allMeanReturns)
            break
        end
        % plot returns histogram
        figure(5)
        h(plotIndex) = subplot(plotRows, plotCols, plotIndex);
        histfit(allReturns(:,plotIndex), 30, 'Normal');
        title(allLabels{plotIndex});
        grid on;
        % plot returns qqplot on another figure
        figure(6)
        % subplot main title
        subplot(plotRows, plotCols, plotIndex);
        normplot(allReturns(:,plotIndex));
        title(allLabels{plotIndex});
        set(gca, 'YTick', [] );
    end
end

%% Plot the sequence of MaxDrawDown for all the assets
figure;
mdd = zeros(size(allReturns));
minMdd = zeros(size(allReturns,2),1);
for i = 1:size(allReturns,2)
    [mddCurr,mddMinCurr] = MDD(allReturns(:,i));
    mdd(:,i) = mddCurr*100;
    minMdd(i) = mddMinCurr*100;
end

% ascending order because mdd is negative
[~, order] = sort(minMdd, 'ascend');
% find EU index position
benchPosition = find(order == size(minMdd,1));

p=plot(dates, mdd(:,order));
datePositions=[-45 -45 -46.2 -45 -46.2 -45];
plotEvents(datePositions, 1);
set(gca, 'XTick', fineDateRange);
h = gca; h.XAxis.TickLabelFormat = 'yyyy';
%xlabel('Date');
ylabel('DrawDown [%]');
title('Sequence of DrawDowns');
grid on;
p(benchPosition).LineWidth = 3;
l=legend(allLabels(order), 'location', 'eastoutside');
title(l,'Min MDD[Incr]','FontWeight','normal');

%% Plot mean vs volatility of each asset
figure;
volatility = sqrt(diag(cov(allReturns)));
scatter(volatility, allMeanReturns, 20, 'm', 'Filled');
% add least squares line
regLine=lsline;
regLine.Color = 'k';
regLine.LineStyle = '--';
hold on
% change color of benchmark
scatter(volatility(end), allMeanReturns(end), 20, 'g', 'Filled');
% plot assets labels
for i = 1:length(allLabels)
    text(volatility(i) + 0.02, allMeanReturns(i), allLabels{i}, 'FontSize', 8);
end
xlabel('Volatility [%]');
ylabel('Expected return [%]');
title('Assets returns vs volatility');
grid on;
clear volatility

%% Plot mean vs MDD for all the assets
figure;
mdd = zeros(length(allLabels),1);
for i = 1:length(allLabels)
    [~,mddMinCurr] = MDD(allReturns(:,i));
    mdd(i) = abs(mddMinCurr);
end
scatter(mdd, allMeanReturns, 20, 'm', 'Filled');
regLine=lsline;
regLine.Color = 'k';
regLine.LineStyle = '--';
hold on;
scatter(mdd(end), allMeanReturns(end), 20, 'g', 'Filled');
for i = 1:length(allLabels)
    text(mdd(i) + 0.3, allMeanReturns(i), allLabels{i}, 'FontSize', 12);
end
hold off;
xlabel('|MaxDD| [%]');
ylabel('Expected return [%]');
title('Assets returns vs |MaxDD|');
grid on;

%% Plot mean vs ES(5%) for all the assets
figure;
es = zeros(length(allLabels),1);
for i = 1:length(allLabels)
    es(i) = ES(allReturns(:,i),0.05);
end
scatter(es, allMeanReturns, 20, 'm', 'Filled');
regLine=lsline;
regLine.Color = 'k';
regLine.LineStyle = '--';
hold on;
scatter(es(end), allMeanReturns(end), 20, 'g', 'Filled');
for i = 1:length(allLabels)
    text(es(i,1) + 0.02, allMeanReturns(i), allLabels{i}, 'FontSize', 12);
end
hold off;
xlabel('|Expected Shortfall| [%]');
ylabel('Expected return [%]');
title('Assets returns vs |Expected Shortfall (5%)|');
grid on;


%% Analysis of the dinamics in the returns series
% ACF of returns
figure;
for i=1:size(allReturns,2)
    subplot(4,5,i);
    % store ACF bounds
    [~,~,bounds]=autocorr(allReturns(:,i));
    % plot ACF values
    autocorr(allReturns(:,i));
    % set y limits for readability (0th lag = 1 will be out of bounds)
    % *3 is just an arbitrarily large value
    axis([0 inf bounds(2)*3 bounds(1)*3])
    ylabel('');
    title(allLabels{i});
end
sgtitle("ACF of assets returns");

% ACF of squared returns, as proxy of variances
figure;
for i=1:size(allReturns,2)
    subplot(4,5,i);
    % store ACF bounds
    [~,~,bounds]=autocorr(allReturns(:,i).^2);
    % plot ACF values
    autocorr(allReturns(:,i).^2);
    % set y limits for readability (0th lag = 1 will be out of bounds)
    % *3 is just an arbitrarily large value
    axis([0 inf bounds(2)*4 bounds(1)*4])
    ylabel('');
    title(allLabels{i});
end
sgtitle("ACF of assets square returns");

% Correlations of assets returns at different lags
j=1;
for lag=[0,1,3,5]
    subplot(2,2,j)
    imagesc(corr_matrix(allReturns,lag));
    caxis([-1 1]); colorbar;
    set(gca,'Ytick',1:size(allReturns,2),'YTickLabel',allLabels,'FontSize',9);
    set(gca,'Xtick',1:size(allReturns,2),'XTickLabel',allLabels,'FontSize',9);
    title(['Assets returns cross correlation [lag ' num2str(lag) ']']);
    j = j+1;
end

% Correlations of assets squared returns at different lags
j=1;
for lag=[0,1,3,5]
    subplot(2,2,j)
    imagesc(corr_matrix(allReturns.^2,lag));
    caxis([-1 1]); colorbar;
    set(gca,'Ytick',1:size(allReturns.^2,2),'YTickLabel',allLabels,'FontSize',9);
    set(gca,'Xtick',1:size(allReturns.^2,2),'XTickLabel',allLabels,'FontSize',9);
    title(['Assets square returns cross correlation [lag ' num2str(lag) ']']);
    j = j+1;
end

%% Summary table with some relevant statistics
Mean=mean(allReturns)';
Median=median(allReturns)';
StDev=sqrt(var(allReturns))';
Min=min(allReturns)';
Max=max(allReturns)';
Skew=skewness(allReturns)';
Kurt=kurtosis(allReturns)';
ES=-es; MaxDD=-mdd;
T=table(Mean,Median,StDev,Min,Max,Skew,Kurt,ES,MaxDD,'RowNames',allLabels);
format bank
T
