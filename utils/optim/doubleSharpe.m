function out=doubleSharpe(w,returns, rollingReturns)
% returns: full returns
% rollingReturns: returns of current rolling window

% bootstrap returns and compute sharpe ratio
% M should be higher
M=100;     % number of simulated paths
N=60;      % length of simulated paths
bs=6;        % average block length

% boostrapping row indices
[indices,~]=stationary_bootstrap_MC((1:size(returns,1))',M,bs,N);

simSharpeRatios=zeros(M,1);
for j=1:M
    % generate simulated returns
    curReturns = returns(indices(:,j),:);
    mR = mean(curReturns);
    vR = cov(curReturns);
    simSharpeRatios(j) = sharpe(w,mR',vR);
end

mR = mean(rollingReturns);
vR = cov(rollingReturns);
sharpeStd = std(simSharpeRatios);

out = sharpe(w,mR',vR)/sharpeStd;
end