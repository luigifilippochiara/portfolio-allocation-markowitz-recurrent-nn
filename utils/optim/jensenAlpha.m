function out=jensenAlpha(w,returns,eqReturns)
portReturns = mean(returns*w);
out=portReturns-eqReturns;
end