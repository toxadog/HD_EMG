clear all
close all
for s=1:10
    for f=30:20:70
        for t=1:2
            clearvars -except s f t

            %Lecture fichier
            fid=fopen(strcat('C:\Users\Mathieu\Desktop\Stage\Données\s',num2str(s),'_e2_f',num2str(f),'_t',num2str(t),'_SplitNum_1.sig'),'r');
            data = fread(fid,[68 inf],'int32','b');

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

            %filtrage
            for j = 50:50:450
                [b a]=butter (2,[j-5 j+5]/(fs/2),'stop');
                for i=1 : 64
                    Emg(:,i) =filtfilt(b,a,Emg(:,i));
                end
            end
            [b a]=butter (8,[20 500]/(fs/2),'bandpass');
            for i=1 : 64
                Emg(:,i) =filtfilt(b,a,Emg(:,i));
            end
            test=zeros(L_data,68);
            test(:,1:64) = Emg;
            test(:,65:68)=mat_Force;
            fw=fopen(strcat('filtre_s',num2str(s),'_e2_f',num2str(f),'_t',num2str(t),'_SplitNum_1.sig'),'w');
            fwrite(fw,test','int32','b');
            
            fclose(fid);
            fclose(fw);
        end
    end
end
