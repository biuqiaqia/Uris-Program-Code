function H_MSL = Get_HK_MSL(Latitude,Longitude,DTM_data) 
% Get HK Mean Sea Level Height by Longitude and Latitude
% Based on HK standard

Longitude_HK80 = Longitude - 8.8/3600;
Latitude_HK80 = Latitude + 5.5/3600;

% Constants
D2R = pi/180;
R2D = 180/pi;
N_0 = 819069.80;
E_0 = 836694.05;
lambda_0 = (114 + 10/60 + 42.80/3600)*D2R;
m_0 = 1;
M_0 = 2468395.723;
a = 6378388;
e2 = 6.722670022e-3;
e4 = e2*e2;
A0 = 1 - (e2/4) - (3*e4)/64;
A2 = 3/8*(e2+e4/4);
A4 = 15/256*e4;

%---HK80 phi%lambda -> HK80 Grid ---
phi = Latitude_HK80*D2R;
lambda = Longitude_HK80*D2R;
nu_s = a / ((1-e2*sin(phi)*sin(phi))^0.5);
rho_s = a*(1-e2)/((1-e2*sin(phi)*sin(phi))^1.5);
psi_s = nu_s/rho_s;

M = a * (A0*phi - A2*sin(2*phi) + A4*sin(4*phi));
N = N_0 + m_0*((M-M_0)+nu_s*sin(phi)*((lambda-lambda_0)*(lambda-lambda_0)/2)*cos(phi));
E = E_0 + m_0*(nu_s*(lambda-lambda_0)*cos(phi)+nu_s*((lambda-lambda_0)^3)/6*(cos(phi))^3*(psi_s-tan(phi)*tan(phi)));

id_E = ceil((E-DTM_data.R.XWorldLimits(1,1))/5) + 1;
id_N = size(DTM_data.Z,1) - ceil((N-DTM_data.R.YWorldLimits(1,1))/5);
H_HKPD = DTM_data.Z(id_N,id_E);
H_MSL = H_HKPD-1.230;






























