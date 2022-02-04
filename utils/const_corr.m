%% return the variance-covariance matrix in a constant correlation model
% input matrix of returns for the assets
function [C] = const_corr(A)
   C = corr(A); n = size(C,1); 
   rho = (sum(C,'all') - n) * (0.5) * (n*(n-1)*0.5);
   sigma = sqrt(var(A));
   for i=1:n
       for j=1:n
           if i==j
               C(i,j) = sigma(i)^2;
           elseif i~=j
               C(i,j) = rho * sigma(i) * sigma(j);
           end;
       end;
   end;