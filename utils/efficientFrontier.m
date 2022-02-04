function [ef] = efficientFrontier(scalars, returns)
% computer scalars used to plot the unconstrained efficient frontier
% Inputs:
% * A,B,C,D scalars computed via mkScalars function
% * input returns of the portfolios
    s = scalars;
    ef = sqrt((s.C/s.D)*(returns.^2)-2*(s.B/s.D)*returns + (s.A/s.D));
end