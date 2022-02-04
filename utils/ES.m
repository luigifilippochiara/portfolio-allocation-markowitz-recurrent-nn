function out=ES(ret,a)
    vara=quantile(ret,a);
    exps=mean(ret(ret<vara));
    out=abs(exps);
end