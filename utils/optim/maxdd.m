function out=maxdd(w,ret)
    % scale returns by 100% because ret is % returns
    ret=ret/100;
    M=mean(ret);
    pret=ret*w;
    dt=zeros(length(pret),1);
    dt(1,1)=min(0,pret(1,1));
    for i=2:length(pret)
       dt(i,1)=min(0,(1+dt(i-1,1))*(1+pret(i))-1);
    end
    ddm=min(dt); % take largest loss
    out=(M*w)/abs(ddm);
end