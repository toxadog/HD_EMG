function Iextrapol = Iextrapolation( I, scale )
% Extrapolation bilin�aire des signaux � chaque instant 
% Input : matrice carr� de taille cot� A et de longueur = nombre
% d'�chantillons temporels
      for i=1:length(I)
          Iextrapol(:,:,i)=imresize(I(:,:,i),scale,'bilinear');
      end

end

