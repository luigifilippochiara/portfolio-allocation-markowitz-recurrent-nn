function [PortWeights, PortReturns]=funGMV(port,returns,j,w)
pwgtGMV=port.estimateFrontierLimits('Min');
PortWeights=pwgtGMV';
PortReturns=returns(j+w,:)*(pwgtGMV);
end