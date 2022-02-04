%% return the variance-covariance matrix in a constant correlation model
% input matrix of returns for the assets
function [C] = const_corr2(A, S)
   B = corr(A); n = size(B,1);
   sigma = sqrt(var(A));
if nargin == 1
   rho = (sum(B,'all') - n) * (0.5) / (n*(n-1)*0.5);
   C = zeros(size(B));
   for i=1:n
       for j=1:n
           if i==j
               C(i,j) = sigma(i).^2;
           elseif i~=j
               C(i,j) = rho * sigma(i) * sigma(j);
           end
       end
   end
else
   k = size(S,2);
   rho = zeros(k,k);
   for i=1:k
       for j=1:k
           if i == j
           rho(i,j) = (sum(B(S(:,i) == 1, S(:,j) == 1),'all') - sum(S(:,i))) * 0.5 / (sum(S(:,i)) * (sum(S(:,i)) - 1) * 0.5); 
           elseif i~=j
           rho(i,j) = sum(B(S(:,i) == 1, S(:,j) == 1),'all') / (sum(S(:,i)) * sum(S(:,j)));
           end
       end
   end
   C = zeros(size(B));
   % create the correlation matrix
   for i=1:k
       for j=1:k
       C(S(:,i) == 1,S(:,j) == 1) = rho(i,j);
       end
   end
   
   % create the covaraince matrix
   for i=1:n
       for j=1:n
           if i==j
               C(i,j) = sigma(i).^2;
           elseif i~=j
               C(i,j) = C(i,j) * sigma(i) * sigma(j);
           end
       end
   end
end
   C = (C+C')*(0.5);
end
   