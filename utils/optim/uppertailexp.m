function out=uppertailexp(w,ret,a)
    pret=ret*w; % port returns
    vara=quantile(pret,a);
    out=mean(pret(pret>=vara));
end