function x=ERCfun(w,Sigma)

% un-standardized Marginal Risks
R = Sigma*w;

% un-standardized Risk Contributions
RC= w.*R;

% squared deviations across all RC pairs
SqRC=(repmat(RC,1,length(RC))-repmat(RC',length(RC),1)).^2;

% computing criterion function
x =sum(sum(SqRC))*0.5;
end