
%% Centre pixels considérés
Emg_propa_centre_1 = (Emg_propa'-repmat(mean(Emg_propa'),length(Emg_propa),1))';
%% Interpolation
clear TF_Emg_propa_centre ZP_Emg_propa_centre X Y
TF_Emg_propa_centre= fftshift(fft (Emg_propa_centre_1'));
ZP_Emg_propa_centre = real(ifft(ifftshift([zeros(20000,3);TF_Emg_propa_centre;zeros(20000,3)]))); % icizeros (len/2,nb colonne)
X = 0:1/fs:(length(Emg_propa_centre_1)-1)/fs;
Y = 0:1/(fs*2):(length(ZP_Emg_propa_centre)-1)/(fs*2);
Emg_propa_centre = ZP_Emg_propa_centre';
%% inter-correlation + mesure retard
clear dist_time
dist_time_2 = zeros(2,length(q_ok));
dist_time_3 = zeros(2,length(q_ok));
for j = 2:length(q_ok)
    clear d_m_corr d_m_corr_pic d_m_corr_pic_2 d_m_corr_pic_3 y x sa sb idx noise_pic std_pic mean_pic
    y=abs(Emg_propa_centre(j-1,:)); % detection pics
    x=1:length(Emg_propa_centre(j-1,:));
    sa=sign(diff([-inf y]));
    sb=sign(diff([-inf fliplr(y)]));
    sb=sb(end:-1:1);
    idx=(sa==1 & sb==1); % possible : & y > quantile(y, 0.80));
    % figure
    % plot(x,y,'b-',x(idx),y(idx),'r*')
    d_m_corr= zeros(2,length(Emg_propa_centre(1,:)));
    for i=81:length(Emg_propa_centre(1,:))-80
        if idx(i) == 1 % donc pic
            [MAX IMAX]=max(xcorr(Emg_propa_centre(j-1,i-80:i+80),Emg_propa_centre(j,i-80:i+80)));
            d_m_corr(1,i)= MAX/std(xcorr(Emg_propa_centre(j-1,i-80:i+80),Emg_propa_centre(j,i-80:i+80))); % a voir
            d_m_corr(2,i) = IMAX-161;
        end
    end
    d_m_corr_pic = d_m_corr;
    d_m_corr_pic (:,idx == 0) = []; %d_m_corr : nombre de points Emg; d_m_corr_pic : slmt_pic
    d_m_corr_pic_2=d_m_corr_pic;
    d_m_corr_pic_2 (:,d_m_corr_pic(1,:)<mean(d_m_corr_pic(1,:))) = []; % seuil = moy du max des intercors 
    
    
    K=2.32;
    noise_pic=zeros(size(d_m_corr_pic_2(1,:)));
    mean_pic=mean(d_m_corr_pic_2(2,:));
    std_pic=std(d_m_corr_pic_2(2,:));
    for k = 1:length(noise_pic)
        if abs(d_m_corr_pic_2(2,k)-mean_pic)>K*std_pic
            noise_pic(k)=1;
        end
    end
    d_m_corr_pic_3 = d_m_corr_pic_2; % sans valeurs abérrantes
    d_m_corr_pic_3 (:,noise_pic==1)= [];
            
    figure;
    subplot(3,2,1);plot(d_m_corr_pic(1,:));title (strcat(num2str(j-1),'-',num2str(j),'max correlation(avant seuillage max intercorr)'))
    subplot(3,2,2);plot(d_m_corr_pic(2,:));title (strcat(num2str(j-1),'-',num2str(j),'retard(avant seuillage max intercorr)'))
    subplot(3,2,3);plot(d_m_corr_pic_2(1,:));title (strcat(num2str(j-1),'-',num2str(j),'max correlation(après seuillage max intercorr)'))
    subplot(3,2,4);plot(d_m_corr_pic_2(2,:));title (strcat(num2str(j-1),'-',num2str(j),'retard(après seuillage max intercorr)'))
    subplot(3,2,5);plot(d_m_corr_pic_3(1,:));title (strcat(num2str(j-1),'-',num2str(j),'max correlation(sans valeurs abérrantes)'))
    subplot(3,2,6);plot(d_m_corr_pic_3(2,:));title (strcat(num2str(j-1),'-',num2str(j),'retard(sans valeurs abérrantes)'))
    
     dist_time_2(2,j-1)= mean(d_m_corr_pic_2(2,:));
     dist_time_3(2,j-1)= mean(d_m_corr_pic_3(2,:));
     dist_time_2(1,j-1) = sqrt((q_ok(1,j-1)-q_ok(1,j))^2+(q_ok(2,j-1)-q_ok(2,j))^2);
     dist_time_3(1,j-1) = sqrt((q_ok(1,j-1)-q_ok(1,j))^2+(q_ok(2,j-1)-q_ok(2,j))^2);
end
%% pour memoire filtre LPD
% for i = 41:30000
% for j = 0:40
% Emg_tt(i)=sin(j*pi/40)*(Emg_propa(1,i+j)-Emg_propa(1,i-j));
% end
% end