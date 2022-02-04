%% funzione per il calcolo di matrici di cross-correlazione ad uno specifico ritardo
function [C] = corr_matrix(A, j)
    % j >= 0...
    C = zeros(size(A,2));
    for i=1:size(C,1)
        for k=1:size(C,1)
            C(i,k) = corr(A((j+1):end,i), A(1:(end-j),k));
        end
    end
