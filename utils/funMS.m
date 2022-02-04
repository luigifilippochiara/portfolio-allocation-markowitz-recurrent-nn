function [PortWeights, PortReturns, pwgtMS]=funMS(p,returns,j,w,pwgtMS,ErS)
testW=p.estimateFrontierLimits('Max');
testRet=ErS(j,:)*testW;
    if testRet>0
        try
            pwgtMS=estimateMaxSharpeRatio(p);
        catch
            disp(['Error while computing MSRatio at iteration ' num2str(j)])
        end
    end
    PortWeights=pwgtMS';
    PortReturns=returns(j+w,:)*(pwgtMS);
end