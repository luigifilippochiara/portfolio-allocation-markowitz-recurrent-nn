function out=sharpeInfoRatio(w,returns,benchReturns)
%R=mean(rI);
%V=cov(rI);

portReturns = returns*w;
excessReturns = portReturns-benchReturns;

out=mean(excessReturns)/std(excessReturns);
end