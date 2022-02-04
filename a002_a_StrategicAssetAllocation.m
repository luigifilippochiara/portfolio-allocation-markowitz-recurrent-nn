%% 2.a Strategic asset allocation
clear
clc
load data.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
    useful plots here /usr/local/MATLAB/R2018b/examples/finance
    portfolioexample.mlx which uses portfolioexamples_plot.m
    portfoliowithturnoverconstraintsexample.m
    turnoverconstraintplot.m
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%}

%% Compute EF over the last 5 years
% restrict to last 5 years

lastDate = dates(end);
firstDate = lastDate - calyears(5);
[datesRg, returnsRg, meanReturnsRg, varcovRg, ~] = restrictRange(firstDate, lastDate, dates, returns);
[~, benchReturnsRg, benchMeanReturnRg, benchVarcovRg, benchRiskRg] = restrictRange(firstDate, lastDate, dates, benchReturns);


%%
C = cbrewer('qual','Paired',20);
%% Portfolio weights, return, risk

nPort = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Sample mean approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init portfolio class
p = Portfolio;
p = p.setAssetList(labels);
% using sample moments
p = p.setAssetMoments(meanReturnsRg, varcovRg);
% no shorting
p = p.setDefaultConstraints;

sm = getplotqnt(nPort,p,meanReturnsRg,varcovRg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Equilibrium approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
showPlots = true;
[meanReturnsEq,varcovEq] = getEquilibrium(returnsRg,meanReturnsRg,benchReturnsRg,benchMeanReturnRg,benchRiskRg,showPlots,labels);

pEq = Portfolio;
pEq = pEq.setAssetList(labels);
pEq = pEq.setAssetMoments(meanReturnsEq, varcovEq);
pEq = pEq.setDefaultConstraints;

eq = getplotqnt(nPort,pEq,meanReturnsEq,varcovEq);

% EW equlibrium returns
lm = regstats(returnsRg*ones(size(returnsRg,2),1)./size(returnsRg,2), benchReturnsRg, 'linear', {'beta', 'r'});
betaEW = lm.beta(2);
% compute equlibrium returns without risk free
EWeqReturns = betaEW*benchMeanReturnRg;
EWeqRisk = sqrt(betaEW^2*benchRiskRg.^2 + var(lm.r));

%% Plot EF over last 5 years
figure;
smPlot = subplot(1,2,1);
portplot(['Sample mean approach' newline datestr(firstDate,'dd/mm/yyyy') ' - ' datestr(lastDate,'dd/mm/yyyy')], ...
    {'line', sm.pRisks, sm.pReturns, {'No-shorting EF'}, '', 1}, ...
    {'line', sm.uncRisks', sm.uncReturns', {'Unconstrained EF'}, '', 1}, ...
    {'scatter', sm.GMVrestr.risk, sm.GMVrestr.return, {'GMV'}}, ...
	{'scatter', sm.TANrestr.risk, sm.TANrestr.return, {'TAN'}}, ...
    {'scatter', sm.GMV.risk, sm.GMV.return, {'GMV'}}, ...
	{'scatter', sm.TAN.risk, sm.TAN.return, {'TAN'}}, ...
    {'scatter', sqrt(diag(varcovRg)), meanReturnsRg, {}, '.k'}, ...
    {'scatter', sm.EWrisk, sm.EWreturn, {'EW'}, 'b'}, ...
	{'scatter', benchRiskRg, benchMeanReturnRg, {'Bench'}, 'b'});
eqPlot = subplot(1,2,2);
portplot(['Equilibrium approach' newline datestr(firstDate,'dd/mm/yyyy') ' - ' datestr(lastDate,'dd/mm/yyyy')], ...
    {'line', eq.pRisks, eq.pReturns, {'No-shorting EF'}, '', 1}, ...
    {'line', eq.uncRisks', eq.uncReturns', {'Unconstrained EF'}, '', 1}, ...
    {'scatter', eq.GMVrestr.risk, eq.GMVrestr.return, {'GMV'}}, ...
	{'scatter', eq.TANrestr.risk, eq.TANrestr.return, {'TAN'}}, ...
    {'scatter', eq.GMV.risk, eq.GMV.return, {'GMV'}}, ...
	{'scatter', eq.TAN.risk, eq.TAN.return, {'TAN'}}, ...
    {'scatter', sqrt(diag(varcovEq)), meanReturnsEq, {}, '.k'}, ...
    {'scatter', EWeqRisk, EWeqReturns, {'EW'}, 'b'}, ...
	{'scatter', benchRiskRg, benchMeanReturnRg, {'Bench'}, 'b'});
%linkaxes([smPlot, eqPlot], 'xy');

%% Weights of TAN portfolio no-shorting vs unrestricted
figure
bar([eq.TAN.weights, eq.TANrestr.weights]);
legend({'Unrestricted', 'Restricted'});
set(gca,'Xtick',1:length(labels),'XTickLabel',labels);
title('Equilibrium TAN portfolio weights');
grid on;

%% Plot blocks of 5 years
% plotting together unconstrained and constrained EF for comparison

lastDate = dates(end);
% get two blocks of 5 years starting from the last date
dateBlocks = NaT(2,2);
for i=1:size(dateBlocks,1)
    dateBlocks(i,2) = lastDate;
    dateBlocks(i,1) = lastDate - calyears(5);
    lastDate = dateBlocks(i,1) - caldays(1);
end

returnsPeriod=zeros(2,60,18);
col = C([floor(linspace(1,size(C,1),size(dateBlocks,1)))],:);
for i=1:size(dateBlocks,1)    
    firstDate = dateBlocks(i,1);
    lastDate = dateBlocks(i,2);
    [datesRg, returnsRg, meanReturnsRg, varcovRg, ~] = restrictRange(firstDate, lastDate, dates, returns);
    [~, benchReturnsRg, benchMeanReturnRg, benchVarcovRg, benchRiskRg] = restrictRange(firstDate, lastDate, dates, benchReturns);
    returnsPeriod(i,:,:)=returnsRg;
    
    % Sample mean approach
    p = Portfolio;
    p = p.setAssetList(labels);
    p = p.setAssetMoments(meanReturnsRg, varcovRg);
    p = p.setDefaultConstraints;
    sm = getplotqnt(nPort,p,meanReturnsRg,varcovRg);

    % Equilibrium approach
    [meanReturnsEq,varcovEq] = getEquilibrium(returnsRg,meanReturnsRg,benchReturnsRg,benchMeanReturnRg,benchRiskRg);
    pEq = Portfolio;
    pEq = pEq.setAssetList(labels);
    pEq = pEq.setAssetMoments(meanReturnsEq, varcovEq);
    pEq = pEq.setDefaultConstraints;
    eq = getplotqnt(nPort,pEq,meanReturnsEq,varcovEq);
    
    if i==1
        figure(5);
        smPlot = subplot(1,2,1);
        portplot('No-shorting SM', ...
            {'line', sm.pRisks, sm.pReturns, {['SM ', datestr(firstDate) ' ' datestr(lastDate)]}, 'r', 1}, ...
            {'scatter', sm.GMVrestr.risk, sm.GMVrestr.return, {'GMV'}, 'r'}, ...
            {'scatter', sm.TANrestr.risk, sm.TANrestr.return, {'TAN'}, 'r'});
        hold on
        eqPlot = subplot(1,2,2);
        portplot('No-shorting EQ', ...
            {'line', eq.pRisks, eq.pReturns, {['EQ ', datestr(firstDate) ' ' datestr(lastDate)]}, 'r', 1}, ...
            {'scatter', eq.GMVrestr.risk, eq.GMVrestr.return, {'GMV'}, 'r'}, ...
            {'scatter', eq.TANrestr.risk, eq.TANrestr.return, {'TAN'}, 'r'});
        hold on
        %linkaxes([smPlot, eqPlot], 'xy');
        figure(6)
        smPlot = subplot(1,2,1);
        portplot('Unconstraned SM', ...
            {'line', sm.uncRisks', sm.uncReturns', {['SM ', datestr(firstDate) ' ' datestr(lastDate)]}, 'r', 1}, ...
            {'scatter', sm.GMV.risk, sm.GMV.return, {'GMV'}, 'r'}, ...
            {'scatter', sm.TAN.risk, sm.TAN.return, {'TAN'}, 'r'});
        hold on
        eqPlot = subplot(1,2,2);
        portplot('Unconstraned EQ', ...
            {'line', eq.uncRisks', eq.uncReturns', {['EQ ', datestr(firstDate) ' ' datestr(lastDate)]}, 'r', 1}, ...
            {'scatter', eq.GMV.risk, eq.GMV.return, {'GMV'}, 'r'}, ...
            {'scatter', eq.TAN.risk, eq.TAN.return, {'TAN'}, 'r'});
        hold on
        %linkaxes([smPlot, eqPlot], 'xy');
    else
        figure(5)
        textDx = 0.005;
        % constrained sm
        subplot(1,2,1)
        plot(sm.pRisks, sm.pReturns, 'Color', col(i,:), 'DisplayName', ['SM ', datestr(firstDate) ' ' datestr(lastDate)])
        scatter(sm.GMVrestr.risk, sm.GMVrestr.return, [], col(i,:), 'Filled', 'HandleVisibility','off')
        text(sm.GMVrestr.risk +textDx, sm.GMVrestr.return, 'GMV', 'Fontsize', 8);
        scatter(sm.TANrestr.risk, sm.TANrestr.return, [], col(i,:), 'Filled', 'HandleVisibility','off')
        text(sm.TANrestr.risk +textDx, sm.TANrestr.return, 'TAN', 'Fontsize', 8);
        % constrained eq
        subplot(1,2,2)
        plot(eq.pRisks, eq.pReturns, 'Color', col(i,:), 'DisplayName', ['EQ ', datestr(firstDate) ' ' datestr(lastDate)])
        scatter(eq.GMVrestr.risk, eq.GMVrestr.return, [], col(i,:), 'Filled', 'HandleVisibility','off')
        text(eq.GMVrestr.risk +textDx, eq.GMVrestr.return, 'GMV', 'Fontsize', 8);
        scatter(eq.TANrestr.risk, eq.TANrestr.return, [], col(i,:), 'Filled', 'HandleVisibility','off')
        text(eq.TANrestr.risk +textDx, eq.TANrestr.return, 'TAN', 'Fontsize', 8);
        
        figure(6)
        % unconstrained sm
        subplot(1,2,1)
        plot(sm.uncRisks', sm.uncReturns', 'Color', col(i,:), 'DisplayName', ['SM ', datestr(firstDate) ' ' datestr(lastDate)])
        scatter(sm.GMV.risk, sm.GMV.return, [], col(i,:), 'Filled', 'HandleVisibility','off')
        text(sm.GMV.risk +textDx, sm.GMV.return, 'GMV', 'Fontsize', 8);
        scatter(sm.TAN.risk, sm.TAN.return, [], col(i,:), 'Filled', 'HandleVisibility','off')
        text(sm.TAN.risk +textDx, sm.TAN.return, 'TAN', 'Fontsize', 8);
        % unconstrained eq
        subplot(1,2,2)
        plot(eq.uncRisks', eq.uncReturns', 'Color', col(i,:), 'DisplayName', ['EQ ', datestr(firstDate) ' ' datestr(lastDate)])
        scatter(eq.GMV.risk, eq.GMV.return, [], col(i,:), 'Filled', 'HandleVisibility','off')
        text(eq.GMV.risk +textDx, eq.GMV.return, 'GMV', 'Fontsize', 8);
        scatter(eq.TAN.risk, eq.TAN.return, [], col(i,:), 'Filled', 'HandleVisibility','off')
        text(eq.TAN.risk +textDx, eq.TAN.return, 'TAN', 'Fontsize', 8);
    end
end

%% Testing the hypothesis that the inputs are the same in both periods

% Log-Ratio Test - last five years inputs equal to the previous period
test = LR(squeeze(returnsPeriod(1,:,:)), mean(squeeze(returnsPeriod(2,:,:))) , cov(squeeze(returnsPeriod(2,:,:))));
1-chi2cdf(test,189) % upper probability
chi2inv(0.95,189) % theoric quantile
% we reject the H0

% compute an alternative quantile(under H0) using bootstrap
% the assumptions of multivariate normal returns and iid are wrong!
[indices,~] = stationary_bootstrap_MC((1:60)',100000,6,60);
test2 = zeros(100000,1);

for j=1:100000
   test2(j,1) =  LR(squeeze(returnsPeriod(2,indices(:,j),:)), mean(squeeze(returnsPeriod(2,:,:))) , cov(squeeze(returnsPeriod(2,:,:))));
end

quantile(test2, 0.95)
% approximate 578.74 against the theoric 222.08
% we reject the hyphotesis that the inputs are stable over time.
% Conclusion: we should use rolling methods to estimate the inputs of our
% strategies. 