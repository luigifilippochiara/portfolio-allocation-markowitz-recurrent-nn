function out=modiglianiRAP(w,returns,eqVarcov)

eqRisks = sqrt(diag(eqVarcov));
portRisk = sqrt(w'*varcov*w);
% leverage param
gammaP = eqRisks/portRisk;

portReturns = returns*w;
out=gammaP*mean(portReturns);
end