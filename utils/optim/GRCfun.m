function x=GRCfun(w,Sigma,b)

% inputs
% w - portfolio weights
% Sigma - covariance of the assets
% b - risk budgets in percentages, to have different contributions to risk

% total risk
R=sqrt((w')*Sigma*w);

% un-standardized Marginal Risks
mR = Sigma*w;

% standardized Risk Contributions
RC= w.*mR/R;

% squared deviations between RC and fractions of total risk
SqRC=(RC - b*R).^2;

% computing criterion function
x =sum(SqRC);
end