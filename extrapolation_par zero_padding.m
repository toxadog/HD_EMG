%              extrapolation par zero-padding
%             TF_8 =  fftshift(fft2(matrice_energie));
%             TF_8(9,2:9)=fliplr(conj(TF_8(1,1:8)));
%             TF_8(2:9,9)=flipud(conj(TF_8(1:8,1)));
%             TF_32 = ifftshift(padarray(TF_8,[12,12]));
%             matrice_energie_32=ifft2 (TF_32);
%             figure (3)
%             imagesc(matrice_energie_32)
%             title ('extrapolation par zero-padding')