function [scalars] = mkScalars(meanReturns, varcov)
    % computes A,B,C,D scalars for standard mk
    % inputs:
    %   r: returns vector in column form
    %   cov: variance covariance matrix
    scalars.A=(meanReturns/varcov)*(meanReturns');
    scalars.B=(meanReturns/varcov)*ones(length(meanReturns),1);
    scalars.C=(ones(1,length(meanReturns))/(varcov))*ones(length(meanReturns),1);
    scalars.D=scalars.A*scalars.C-scalars.B*scalars.B;
end