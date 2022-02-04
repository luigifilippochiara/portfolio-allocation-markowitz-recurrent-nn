function [b]=beta(returns,benchReturns)
    X=[ones(size(returns,1),1) benchReturns];
    b=(X'*X)\(X'*returns);
    b=b(2);
end