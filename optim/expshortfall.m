function out=expshortfall(w,ret,a)
    % M=mean(ret);
    pret=ret*w; % port returns
    vara=quantile(pret,a); % a% quantile of the port return (a value at risk)
    exps=mean(pret(pret<vara)); % this is the exp shortfall
    out=median(pret)/abs(exps);
end