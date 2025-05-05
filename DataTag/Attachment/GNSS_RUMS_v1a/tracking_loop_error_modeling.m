function [sigma_DLL,sigma_PLL,sigma_FLL] = tracking_loop_error_modeling(mea_CN0,lambda,para_cl)

% ------------------------------------------------------------------------
% GNSS DLL,PLL,FLL tracking loop introduced error
%
% based on 'Kaplan, E., & Hegarty, C. (2005). Understanding GPS: principles...
%           and applications. Artech house.'
%
% by GH.Zhang 2020/08/01
% guo-hao.zhang@connect.polyu.hk
% ------------------------------------------------------------------------
% clc;
% clear;
% close all;

goGNSS;
goGNSS.V_LIGHT;

b = 2; %normalized bandwidth
D = 1; %E-L correlator spacing
B_n_DLL = 0.2; %DLL noise bandwidth
T = 0.02; %Predetection integration time
R_c = 1.023e6; %C/A L1 chip rate
L_c = goGNSS.V_LIGHT/R_c; %chip-length
T_c = 1/R_c; %chip period
B_fe = b/T_c; %front end bandwidth

B_n_PLL = 15; %Carrier loop noise bandwidth
% f_L = 1575.42e6;
% lambda = goGNSS.V_LIGHT/f_L;
f_L = goGNSS.V_LIGHT/lambda;%L-band frequency
AD = 1e-10;%Allan deviation (degree)
LOS_jds_max = 98; %maximum LOS jerk dynamic stress (m/s^3)

% id = 1;
% for CN0 = 15:0.1:45
CNR = 10^(mea_CN0/10);

%% DLL
% DLL thermal noise for BPSK-R
if D>=1/b*pi
    sigma_tDLL = sqrt(B_n_DLL/2/CNR*D*(1+2/(T*CNR*(2-D))));
elseif D>1/b
    sigma_tDLL = sqrt(B_n_DLL/2/CNR*(1/B_fe/T_c+B_fe*T_c/(pi-1)*(D-1/B_fe/T_c)^2)*(1+2/(T*CNR*(2-D))));
else
    sigma_tDLL = sqrt(B_n_DLL/2/CNR*(1/B_fe/T_c)*(1+1/T/CNR));
end
% DLL dynamic stress (unaided 3rd order C/A code DLL)
R_e = 0; % neglected by carrier-aided code technique
% DLL noise (in meter)
sigma_DLL = (sigma_tDLL + R_e/3)*L_c;

%% PLL
% PLL thermal noise
sigma_tPLL = 360/2/pi*sqrt(B_n_PLL/CNR*(1+1/(2*T*CNR)));
% Vibration-induced oscillator phase noise
sigma_v = 90.265*sqrt(0.005*(1/20-1/2000));
% Allan Deviation Oscillator Phase Noise (3rd order PLL)
cita_A = 160*AD*f_L/B_n_PLL;
% Dynamic Stress Error
d3R_over_dt3 = LOS_jds_max*360*f_L/goGNSS.V_LIGHT;
cita_e = 0.4828*d3R_over_dt3/B_n_PLL^3;
% PLL_noise (in meter)
sigma_PLL = (sqrt(sigma_tPLL^2+sigma_v^2+cita_A^2)+cita_e/3)/360*lambda;

%% FLL
B_n_FLL = B_n_PLL;
% FLL thermal noise
threshold = 1/12/T;
if CNR >= threshold*10
    F = 1;
else
    F = 2;
end
sigma_tFLL = 1/2/pi/T*sqrt(4*F*B_n_FLL/CNR*(1+1/(T*CNR)));
% Dynamic Stress Error (2nd order FLL)
omega_0 = B_n_FLL/0.53;
f_e = 1/360/(omega_0^2)*d3R_over_dt3;
% FLL_noise (in Hz)
sigma_FLL = sigma_tFLL + f_e/3;

% Noise Threshold
sigma_DLL = min([sigma_DLL,para_cl/6]);
sigma_PLL = min([sigma_PLL,15/360*lambda]);
sigma_FLL = min([sigma_FLL,threshold]);

%% check error
% sigma(id,:) = [CN0,sigma_DLL,sigma_PLL,sigma_FLL];
% id=id+1;
% end

% figure(1);
% grid on;
% hold on;
% plot(sigma(:,1),sigma(:,2),'k-','LineWidth',2);
% xlabel('C/N_0 (dB-Hz)');
% ylabel('DLL 1-\sigma noise (m)');
% ylim([0,25]);
% xlim([15,45]);
% 
% figure(2);
% grid on;
% hold on;
% plot(sigma(:,1),sigma(:,3),'k-','LineWidth',2);
% xlabel('C/N_0 (dB-Hz)');
% ylabel('PLL 1-\sigma noise (m)');
% ylim([0,0.05]);
% xlim([15,45]);
% 
% figure(3);
% grid on;
% hold on;
% plot(sigma(:,1),sigma(:,4),'k-','LineWidth',2);
% xlabel('C/N_0 (dB-Hz)');
% ylabel('FLL 1-\sigma noise (Hz)');
% ylim([0,25]);
% xlim([15,45]);

















