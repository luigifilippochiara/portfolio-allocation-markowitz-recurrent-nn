function [dt, out] =MDD(ret)
    ret=ret/100;
    dt=zeros(length(ret),1);
    dt(1,1)=min(0,ret(1,1));
    for i=2:length(ret)
       dt(i,1)=min(0,(1+dt(i-1,1))*(1+ret(i,1))-1);
    end
    out=min(dt)*100;
end