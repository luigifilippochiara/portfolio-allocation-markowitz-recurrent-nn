function [returns, risk] = estimatePort(meanReturns, varcov, weights)
% computes returns and risk of a portfolio with given weights
% Inputs:
% * mean returns of the assets
% * variance covariance matrix of the assets
% * custom weights
    returns=sum(meanReturns*weights);
    risk=sqrt((weights')*varcov*weights);
end

