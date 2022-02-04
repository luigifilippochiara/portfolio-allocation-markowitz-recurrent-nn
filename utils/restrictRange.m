function [datesRg, returnsRg, meanReturnsRg, varcovRg, riskRg] = restrictRange(firstDate, lastDate, dates, returns)
% recomputes mean and covariances for restricted periods of time
% Inputs;
% * first and last date in datetime format
% * dates array in datetime format
% * asset returns matrix
    dataFilter = (dates >= firstDate) & (dates <= lastDate);
    datesRg = dates(dataFilter);
    returnsRg = returns(dataFilter, :);
    meanReturnsRg = mean(returnsRg);
    varcovRg = cov(returnsRg);
    riskRg = sqrt(diag(varcovRg));
end