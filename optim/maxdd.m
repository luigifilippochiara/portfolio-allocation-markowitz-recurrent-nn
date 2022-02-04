function out=maxdd(w,ret)
    %M=mean(ret);
    ret=ret./100;
    pret=ret*w;
    dt=zeros(length(pret),1);
    dt(1,1)=min(0,pret(1,1));
    for i=2:length(pret)
       dt(i,1)=min(0,(1+dt(i-1,1))*(1+pret(i))-1);
    end
    ddm=min(dt); % take largest loss
    out=median(pret)/abs(ddm);
end