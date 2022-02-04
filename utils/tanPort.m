function [TAN] = tanPort(meanReturns, varcov)
% computes Tangency or Maximum Tradeoff portfolio
% Inputs:
% * mean returns of the assets
% * variance covariance matrix of the assets
    TAN.weights=((varcov)\meanReturns')/sum((varcov)\meanReturns');
    TAN.return=sum(meanReturns*TAN.weights);
    TAN.risk=sqrt((TAN.weights')*varcov*TAN.weights);
end