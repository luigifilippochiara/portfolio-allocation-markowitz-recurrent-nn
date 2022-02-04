function out=rewardToVAR(w,returns)

%varcov = cov(returns);
mPortReturn = mean(returns*w);
%portRisk = sqrt(w'*varcov*w);
% Value-at-Risk
RiskThreshold=0.05;
VAR=mPortReturn./abs(quantile(returns,RiskThreshold));
% should be equal to
% VAR = portvrisk(mPortReturn,portRisk,RiskThreshold);

out=mPortReturn/VAR;
end