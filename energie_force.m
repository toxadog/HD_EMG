clear all
close all
m_rms = [0 0 0]
for S=1:10
    for f=30:20:70
        for t=1:2
            clearvars -except S f t m_rms
%             S=1;
%             f=50;
%             t=2;
            longueur_ech_active=50000;
            
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
                        if i == 1 && j == 1
                            Emg_3d(i,j,:) = 0.5 * (Emg_3d(i+1,j,:) + Emg_3d(i,j+1,:));
                        elseif i == 1 && j == 8
                            Emg_3d(i,j,:) = 0.5 * (Emg_3d(i+1,j,:) + Emg_3d(i,j-1,:));
                        elseif i == 8 && j == 1
                            Emg_3d(i,j,:) = 0.5 * (Emg_3d(i-1,j,:) + Emg_3d(i,j+1,:));
                        elseif i == 8 && j == 8
                            Emg_3d(i,j,:) = 0.5 * (Emg_3d(i-1,j,:) + Emg_3d(i,j-1,:));
                        elseif i==1
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
%% moyenne rms
            if (f == 30)
                m_rms(1) = m_rms(1) + mean(mean(matrice_energie));
            end
            if (f == 50)
                m_rms(2) = m_rms(2) + mean(mean(matrice_energie));
            end
            if (f == 70)
                m_rms(3) = m_rms(3) + mean(mean(matrice_energie));
        
            end
            
        end
    end
end
%% pourcentage correspondant
m_p(1) = m_rms(1)*0.7/m_rms(3);
m_p(2) = m_rms(2)*0.7/m_rms(3);
m_p(3) =0.7;

