function [turnover]=approxTurnover(w)
    w = squeeze(w);
    % subtract each row with previous one
    wDelta = w-lagmatrix(w,1);
    % remove nan rows
    filter = any(isnan(wDelta),2);
    wDelta = wDelta(filter==0,:);
    % sum abs deltas (columns) and divide by 2
    turnover = sum(abs(wDelta), 2)*0.5;
end