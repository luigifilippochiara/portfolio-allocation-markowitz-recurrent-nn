%% funzione che dati i pesi ed una matrice per i gruppi restituisca i pesi aggregati per gruppo
function [W] = w_groups(w, G)
   W = zeros(size(w,1),size(G,2));
   for j=1:size(G,2)
      W(:,j) = sum( w(:,G(:,j) == 1) ,2); 
   end
end