function out=sortino(w,ret)
    M=mean(ret);
    pret=ret*w; % historical return of the portfolio returns with weights w
    pret=pret(pret<0); % keep negative port returns
    out=(M*w)/sqrt(var(pret)); % port returns divided by downside deviations
end