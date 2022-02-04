function [returns] = computeReturns(data)
    % computes returns as % differences between subsequent periods
    returns=((data(2:end,:)./data(1:end-1,:))-1)*100;
end