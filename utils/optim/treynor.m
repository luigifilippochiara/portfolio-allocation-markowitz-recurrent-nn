function out=treynor(w,R,Rb)
    % compute beta on the European market index
    returns=R;
    benchReturn=Rb;
    meanBenchReturn=mean(benchReturn);
    portReturn=returns*w;
    meanPortReturn=mean(portReturn);
    portBeta=(((benchReturn-meanBenchReturn)')*(portReturn-meanPortReturn))...
              /((benchReturn-meanBenchReturn)'*(benchReturn-meanBenchReturn));
    out=meanPortReturn./portBeta;
end