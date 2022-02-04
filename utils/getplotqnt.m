function [r] = getplotqnt(nPort,p,meanReturns,varcov)

% GMV and TAN of constrained portfolio
GMVrestr.weights = p.estimateFrontierLimits('Min');
[GMVrestr.risk, GMVrestr.return] = p.estimatePortMoments(GMVrestr.weights);
TANrestr.weights = p.estimateMaxSharpeRatio();
[TANrestr.risk, TANrestr.return] = p.estimatePortMoments(TANrestr.weights);

% EF of constrained portfolio
pWeights = p.estimateFrontier(nPort);
[pRisks,pReturns] = p.estimatePortMoments(pWeights);

% GMV and TAN of unrestricted portfolio
% The portfolio class optimization algorithm requires portfolio set
% to have at least a finite lower bound. Using handmande approach since we
% are working with unconstrained EF
GMV = gmvPort(meanReturns, varcov);
TAN = tanPort(meanReturns, varcov);
scalars = mkScalars(meanReturns, varcov);

% unconstrained EF
returnLimit = TAN.return*2;
uncReturns = linspace(GMV.return, returnLimit, nPort);
uncRisks = efficientFrontier(scalars, uncReturns);

% equally weighted portfolio
EWweights = ones(1, length(meanReturns))/length(meanReturns);
[EWreturn, EWrisk] = estimatePort(meanReturns, varcov, EWweights');

% storing relevant quantities
r.pRisks = pRisks;
r.pReturns = pReturns;
r.uncRisks = uncRisks;
r.uncReturns = uncReturns;
r.EWrisk = EWrisk;
r.EWreturn = EWreturn;
r.GMV = GMV;
r.TAN = TAN;
r.GMVrestr = GMVrestr;
r.TANrestr = TANrestr;

end

