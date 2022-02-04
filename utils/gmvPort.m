function [GMV] = gmvPort(meanReturns, varcov)
% computes Global Minimum Variance portfolio
% Inputs:
% * mean returns of the assets
% * variance covariance matrix of the assets
    GMV.weights=((varcov)\ones(size(varcov,1),1))/sum((varcov)\ones(size(varcov,1),1));
    GMV.return=sum(meanReturns*GMV.weights);
    GMV.risk=sqrt((GMV.weights')*varcov*GMV.weights);
end

