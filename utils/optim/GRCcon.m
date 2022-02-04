function [c,ceq] = GRCcon(w)
global Sigma;
global Gmat;
global bvec;
% inputs
% w: k-by-1 vector of asset weigths
% Sigma: covariance
% G: q-by-k matrix of q groups
% b: q-by-1 risk budgets for each group
% limiting case G=eye(k) is a risk budget for each asset

% non-linear inequalities are not present so fix them to satisfy the
% fmincon requirement
c=[];

% total risk
R=sqrt((w')*Sigma*w);

% un-standardized Marginal Risks
mR = Sigma*w; % 18x1

% standardized Risk Contributions
RC= w.*mR/R; % 18x1
RC=RC./sum(RC);

% Compute nonlinear equalities
ceq = (Gmat*RC-bvec'*R).^2;   
    
end