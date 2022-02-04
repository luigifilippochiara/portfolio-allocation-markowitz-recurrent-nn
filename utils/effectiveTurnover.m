function [turnover]=effectiveTurnover(w,returns)
    wEnd = zeros(size(w,1)-1, size(w,2));
    W = 1;
    wBegin = w(1,:);
    returns = returns/100;
    for i=1:size(returns,1)-1
        W = W*(1+w(i,:)*returns(i,:)');
        %W = vBegin * (1+returns(i,:))';
        wEnd(i,:) = wBegin.*(1+returns(i,:))/W;
        wBegin = W*w(i+1,:);
    end
    turnover = sum(abs(w(2:end,:)-wEnd),2).*0.5;
end