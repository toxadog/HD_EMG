 %% test watershed : marche pas bien
% I = mat2gray(matrice_energie_32);
% I = histeq(I);
% hy = fspecial('sobel');
% hx = hy';
% Ix = imfilter(I,hy,'replicate');
% Ix = imfilter(I,hx,'replicate');
% Iy = imfilter(I,hy,'replicate');
% g = sqrt(Ix.^2+Iy.^2);
% se = ones(3,3);
% Io = imopen(g,se);
% Ic = imclose(Io,se);
% 
% test = watershed(Ic);
% figure
% imagesc(test)
%% test h-dome : donne les mêmes résultats que moi
% marker = I - quantile(reshape(I,1,32*32), 0.10);
% reconstr = imreconstruct(marker,I);
% test2 = I-reconstr;
% thresh = graythresh(test2);
% test2 = test2 >thresh;
% figure
% imagesc(test2)