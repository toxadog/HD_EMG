%% A ajouter/modifier
% cas defectueux dans un coin : voir energie_force.m
% verif dans calcul seuil : i > 50
% filtrage
% watershed
%% 



clear all
close all
% for S=1:10
%     for f=30:20:70
%         for t=1:2
%             clearvars -except S f t
            S=1;
            f=50;
            t=2;
             longueur_ech_active=40000;
%%            
            %Lecture fichier
            fid=fopen(strcat('C:\Users\Mathieu\Desktop\Stage\Edm_filtre\filtre_s',num2str(S),'_e2_f',num2str(f),'_t',num2str(t),'_SplitNum_1.sig'),'r');
            data = fread(fid,[68 inf],'int32','b');
%%
            %Paramètres
            L_data = length(data);
            sampling_period = 409.6E-6;
            fs = 2441.4;
            v_temps=(0 : sampling_period : (L_data-1)*sampling_period)';

            %Données de travail
            mat_Force = data (65:68,:)';
            Force_index = mat_Force(:,1);
            Force_edm = mat_Force(:,4);
            Force=Force_edm;
            Emg = data (1:64,:)';
%% sert plus
            %filtrage
%             for j = 50:50:450
%                 [b a]=butter (2,[j-5 j+5]/(fs/2),'stop');
%                 for i=1 : 64
%                     Emg(:,i) =filtfilt(b,a,Emg(:,i));
%                 end
%             end
%             [b a]=butter (8,[20 500]/(fs/2),'bandpass');
%             for i=1 : 64
%                 Emg(:,i) =filtfilt(b,a,Emg(:,i));
%             end
%%
            %Détection du mouvement
            % Force = Force - mat_Force(:,3);
            seuil = mean(abs(Force));
            Test = abs(Force)> seuil;
            for i=1:L_data; %suppression des valeurs de forces dûes à la recalibration
                if Test (i) ==1
                    if Test(i-50) ~=1 && Test(i+50) ~=1 
                        Test (i)=0;
                    end
                end
            end
%             essai = find(Test,longueur_ech_active,'first')
%             Emg_active = Emg(essai(1):essai(longueur_ech_active),:);
             Emg_active = Emg;
             Emg_active(Test==0,:)=[];
             Emg_active (longueur_ech_active+1:end,:)=[];
%%
            % Représentation force et seuil
            figure (1)
            subplot (3,1,1)
            title('seuillage de la force')
            se=ones(1,L_data)*seuil;
            plot(abs(Force));
            hold
            plot(se)
            plot(Test*1e+5)
            subplot (3,1,2)
            title('Emg totales')
            plot(Emg(:,10:11));
            subplot (3,1,3)
            title('Emg actives')
            plot(Emg_active(:,10:11));          
%%
            %Detection pixel defectueux
            Emg_non_active = Emg;
            Emg_non_active(Test==1,:)=[];
            
            v_noise = rms(Emg_non_active);
            v_noise(47)=v_noise(46); % pixel défectueux déjà connu
            
            K=2.32;
            noise_en=zeros(size(v_noise));
            mean_en=mean(v_noise);
            std_en=std(v_noise);
            for i = 1:length(v_noise)
                if abs(v_noise(i)-mean_en)>K*std_en
                    noise_en(i)=1;
                end
            end
            
            for i=1:8
                if mod(i,2 )==0
                    m_noise (i,:)= noise_en(64-8*i+1 : 64-8*(i-1));
                else
                    m_noise (i,:)= fliplr(noise_en(64-8*i+1 : 64-8*(i-1)));
                end
            end
%%
           % Création matrice 3D
           for i=1:8
                if mod(i,2 )==0
                   Emg_3d(i,:,:)= Emg_active(:,64-8*i+1 : 64-8*(i-1))';
                else
                   Emg_3d (i,:,:)= fliplr(Emg_active(:,64-8*i+1 : 64-8*(i-1)))';
                end
           end
%%
            %Correction pixel defectueux
%           Correction capteur défectueux
            Emg_3d(3,2,:) = 0.25 * (Emg_3d(4,2,:)+Emg_3d(2,2,:) + Emg_3d(3,1,:)+Emg_3d(3,3,:)); 
            for i=1:8
                for j=1:8
                    if m_noise(i,j)==1 
                        if i==1
                            Emg_3d(i,j,:) = 0.33 * (Emg_3d(i+1,j,:) + Emg_3d(i,j-1,:)+Emg_3d(i,j+1,:));
                        elseif i==8
                            Emg_3d(i,j,:) = 0.33 * (Emg_3d(i-1,j,:) + Emg_3d(i,j-1,:)+Emg_3d(i,j+1,:));
                        elseif j==1
                            Emg_3d(i,j,:) = 0.33 * (Emg_3d(i+1,j,:)+Emg_3d(i-1,j,:)+Emg_3d(i,j+1,:));
                        elseif j==8
                            Emg_3d(i,j,:) = 0.33 * (Emg_3d(i+1,j,:)+Emg_3d(i-1,j,:)+Emg_3d(i,j-1,:));
                        else
                            Emg_3d(i,j,:) = 0.25 * (Emg_3d(i+1,j,:)+Emg_3d(i-1,j,:) + Emg_3d(i,j-1,:)+Emg_3d(i,j+1,:));
                        end
                    end
                end
            end
%%
            % Creation test 
            for i=1:8
                for j=1:8
                    matrice_energie (i,j)= rms (Emg_3d (i,j,:)); 
                end
            end
%             figure (2)
%             imagesc(matrice_energie)
%             title('energie')            
%% sert plus
% %             % Creation matrice energie 
% %             v_energie = rms(Emg_active);
% %             v_energie(47)=v_energie(46);
% %             for i=1:8
% %                 if mod(i,2 )==0
% %                     matrice_energie (i,:)= v_energie(64-8*i+1 : 64-8*(i-1));
% %                 else
% %                     matrice_energie (i,:)= fliplr(v_energie(64-8*i+1 : 64-8*(i-1)));
% %                 end
% %             end
% %             figure (3)
% %             imagesc(matrice_energie)
% %             title('energie')
            
%%
            % extrapolation par zero-padding
%             TF_8 =  fftshift(fft2(matrice_energie));
%             TF_8(9,2:9)=fliplr(conj(TF_8(1,1:8)));
%             TF_8(2:9,9)=flipud(conj(TF_8(1:8,1)));
%             TF_32 = ifftshift(padarray(TF_8,[12,12]));
%             matrice_energie_32=ifft2 (TF_32);
%             figure (3)
%             imagesc(matrice_energie_32)
%             title ('extrapolation par zero-padding')
            % extrapolation par bilinéaire
%             figure (4)
              matrice_energie_32 = imresize(matrice_energie,4,'bilinear');
%             imagesc(matrice_energie_32)
%             title ('extrapolation bilinéaire')
%%          
%             test reconnaissance direction
              
              seuil = quantile(reshape(matrice_energie_32,1,32*32), 0.90);
              logical = matrice_energie_32<seuil;
              logical = 1-logical;
%               possible 
%               se = ones(3,3);
%               logical = imopen(logical,se);
%               logical = imclose(logical,se);
%               figure 
%               imagesc(logical)
              s = regionprops(logical, 'Orientation', 'MajorAxisLength','MinorAxisLength', 'Eccentricity','Centroid');
              phi = linspace(0,2*pi,50);
              cosphi = cos(phi);
              sinphi = sin(phi);
              xbar = s.Centroid(1);
              ybar = s.Centroid(2);
              a = s.MajorAxisLength/2;
              b = s.MinorAxisLength/2;
              theta = pi*s.Orientation/180;
              R = [ cos(theta)   sin(theta); -sin(theta)   cos(theta)];
              xy = [a*cosphi; b*sinphi];
              xy = R*xy;
              x = xy(1,:) + xbar;
              y = xy(2,:) + ybar;
              %test
              %4*cos(-s.Orientation*pi/180)
                   q = 0 :5.66*cos(-s.Orientation*pi/180):s.MajorAxisLength*cos(-s.Orientation*pi/180);
                   q(2,:)= q(1,:)*tan(-s.Orientation*pi/180);
                   q(1,:) = q(1,:)+s.Centroid(1)-median(q(1,:));
                   q(2,:) = q(2,:)+s.Centroid(2)-median(q(2,:));
%               %
%               figure 
%               imagesc(matrice_energie_32)
%               hold on
%               plot(x,y,'r','LineWidth',2);
%               plot(q(1,:),q(2,:),'r','LineWidth',2);
%               print(strcat('Direction_s',num2str(S),'_e2_f',num2str(f),'_t',num2str(t)),'-djpeg','-r0')

%%  Interpolation signaux
    
      for i=1:length(Emg_active)
          Emg_3d_32(:,:,i)=imresize(Emg_3d (:,:,i),4,'bilinear');
      end
%% 
%       Matrice propagation
        q_ok=q;
        [row col]=find(q>32); %ajout un peu sale (il arrive que notre droite de propa sorte de l'image)
        q_ok(row,col) = floor(q_ok(row,col));
        [row col]=find(q_ok>32);
        q_ok(:,col) = [];
        [row col]=find(q<1);
        q_ok(:,col) = [];
        Emg_propa =zeros(length(q_ok),length(Emg_active));
        for i=1:length(q_ok)
            Emg_propa(i,:)=squeeze(Emg_3d_32(int8(q_ok(2,i)),int8(q_ok(1,i)),:));
        end
        figure(8)
        plot(Emg_propa(1,:)');
        hold on
        for i = 2:length(q_ok)
            plot(Emg_propa(i,:)'+i*12e+4);
            drawnow
        end
        title('Emg direction muscle')
  %% 
%       Matrice propagation Emg_diff
        
        for i=1:length(q_ok)-1
            Emg_propa_diff(i,:)=Emg_propa(i,:)-Emg_propa(i+1,:);
        end
        figure(9)
        plot(Emg_propa_diff(1,:)');
        hold on
        for i = 2:length(q_ok)-1
            plot(Emg_propa_diff(i,:)'+i*4e+4);
            drawnow
        end    
        title('Emg diff direction muscle')
%% 
%      Affichage pixel considérés
       matrice_energie_test=matrice_energie_32;
       for i=1:length(q_ok)
           matrice_energie_test(int8(q_ok(2,i)),int8(q_ok(1,i)))=0;
       end
       figure (12)
       imagesc(matrice_energie_test);
       print(strcat('pixel_s',num2str(S),'_e2_f',num2str(f),'_t',num2str(t)),'-djpeg','-r0')
%%
            fclose(fid)
%           close all

%         end
%     end
% end
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