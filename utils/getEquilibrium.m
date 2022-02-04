function [meanReturnsEq,varcovEq] = getEquilibrium(returns,meanReturns,benchReturns,benchMeanReturn,benchRisk,showplots,labels)

% set default values for last 2 variables
if nargin < 6
    showplots = false;
    labels = {};
end

% CAPM estimation

% initialize relevant quantities of CAPM model
assetN = size(returns, 2);
alpha = zeros(assetN, 1);
beta = zeros(assetN, 1);
r2 = zeros(assetN, 1);
pval = zeros(2, assetN);
residuals = zeros(size(returns));
eqReturns = zeros(assetN, 1);

% market equilibrium return and variance
marketEqReturn = benchMeanReturn;
marketVar = benchRisk.^2;

% estimate parameters for each asset
for i = 1:assetN
    % linear model to estimate alpha and beta params
    lm = regstats(returns(:,i), benchReturns, 'linear', {'beta','r','rsquare','tstat'});
    alpha(i)       = lm.beta(1);
    beta(i)        = lm.beta(2);
    r2(i)          = lm.rsquare;
    pval(:,i)      = lm.tstat.pval;
    residuals(:,i) = lm.r;
    % compute equlibrium returns without risk free
    eqReturns(i)   = beta(i)*marketEqReturn;
end

% mean returns and varcov matrix for equilibrium approach
meanReturnsEq=eqReturns';
varcovEq=beta*(beta')*marketVar+diag(diag(cov(residuals)));

% fit plots
if(showplots)
    % plot significance level of parameters and r-squared
    figure
    subplot(1,2,1);
    bar(r2);
    title('R-squared');
    set(gca,'Xtick',1:assetN,'XTickLabel',labels);
    grid on
    subplot(1,2,2);
    bar(pval');
    hold on
    hline = refline(0, 0.05);
    hline.LineStyle = '--';
    hline.Color = 'k';
    set(gca,'Xtick',1:assetN,'XTickLabel',labels);
    title('P-values');
    grid on
    legend({'alpha','beta'});
    % does not work
    %sgtitle('Linear regression fit analysis');

    % plot beta and mean return values
    figure
    subplot(1,2,1);
    bar(beta);
    title('Beta parameters');
    set(gca,'Xtick',1:assetN,'XTickLabel',labels);
    grid on
    subplot(1,2,2);
    bar([meanReturns' eqReturns]);
    ylabel('Return [%]');
    set(gca,'Xtick',1:assetN,'XTickLabel',labels);
    title('Asset returns');
    legend({'Returns with sample mean', 'Returns with CAPM'});
    grid on;
    % does not work
    %sgtitle('Mean returns comparison');
end


end

