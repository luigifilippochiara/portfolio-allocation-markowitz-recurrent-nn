function out=sharpe(w,R,V)
out=(w'*R)/sqrt(w'*V*w);
end