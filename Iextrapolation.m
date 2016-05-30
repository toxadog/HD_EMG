function Iextrapol = Iextrapolation( I, scale )
% Extrapolation bilinéaire des signaux à chaque instant 
% Input : matrice carré de taille coté A et de longueur = nombre
% d'échantillons temporels
      for i=1:length(I)
          Iextrapol(:,:,i)=imresize(I(:,:,i),scale,'bilinear');
      end

end

