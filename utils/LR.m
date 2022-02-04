%% Likelihood-Ratio Test
function [test] = LR(Y, MU, SIGMA)
   r = size(Y,1); SIGMA = SIGMA .* ( (r-1) / r );
   MM = mean(Y);
   VV = cov(Y) .* ( (r-1) / r );
   test = 2 * (sum(log(mvnpdf(Y,MM,VV))) - sum(log(mvnpdf(Y,MU,SIGMA)))); 
end