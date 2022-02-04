function [d]=diversification(w,c)
    w=squeeze(w);
    % weight of EW port
    wEW = (1/c)*ones(size(w));
    % compute diversification wrt EW portfolio
    d=sum(abs(w-wEW), 2);
end