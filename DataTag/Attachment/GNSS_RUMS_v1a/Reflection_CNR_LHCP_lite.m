function CNR_r = Reflection_CNR_LHCP_lite(R_LHCP,CNR0,range,delay)

% calculating LHCP signal CNR
% Referenced by Taro suggested paper
%
% Jia, Yan, and Yuekun Pei. "Remote Sensing in Land Applications by Using 
% GNSS-Reflectometry." Recent Advances and Applications in Remote Sensing. 
% IntechOpen, 2018.
%
% Angelo Joseph. "Measuring_GNSS_Signal_Strength," Inside GNSS 

% [~,idr] = min(refl_dist);%selecting shortest path
% 
% c = norm(SV_pos-user_pos);
% a = norm(SV_pos-refl_pos(idr,:));
% b = norm(user_pos-refl_pos(idr,:));
% 
% % BW = 10*log10(RF_BW);%bandwidth of observation (dB)
% % SNR0 = CNR0 - BW;%open-sky snr in dB
% % SNR0_ = 10^(SNR0/10); %(dBW)
% % 
% % SNR_r_ = SNR0_*(c/(a+b))^2*R_LHCP(idr)^2;
% % CNR_r = 10*log10(SNR_r_) + BW;
% 
% CNR_r = CNR0 + 10*log10((c/(a+b)*abs(R_LHCP(idr)))^2);
% 
% dist_r = refl_dist(idr);
% pos_r = refl_pos(idr,:);

CNR_r = CNR0 + 10*log10((range/(range+delay)*abs(R_LHCP))^2);

