function sigma_cnr = CNR_variance_modeling(CNR)

% Modeling CNR noise with Narrowband-Wideband method
% using curve fitting miu_p & sigma(P_nw)
% input CNR in unit of dB-Hz
% output sigma_cnr in unit of Hz
%
% 'Groves, P. D. (2005). GPS Signal?to?Noise Measurement in Weak Signal and 
%  High?Interference Environments. Navigation, 52(2), 83-94.'

M = 5;
K = 10;
T = 0.02;

cnr = 10^(CNR/10);% in Hz
miu_p = (5.572e-14*cnr^2 + 5*cnr + 50) / (cnr + 50);
sigma_P_nw = 0.4206*exp((-0.002776)*cnr) + 0.1004*exp((-0.0003015)*cnr);
sigma_cnr = (1/T)*(M-1)/((M-miu_p)^2)*sigma_P_nw/sqrt(K);





