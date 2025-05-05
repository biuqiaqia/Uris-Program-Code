function pos = NED_to_ECEF_pos(L_b,lambda_b,h_b)

% Parameters
R_0 = 6378137; %WGS84 Equatorial radius in meters
e = 0.0818191908425; %WGS84 eccentricity
D2R = pi/180;
% Begins
lambda_b=lambda_b*D2R;
L_b=L_b*D2R;
% Calculate transverse radius of curvature using (2.105)
R_E = R_0 / sqrt(1 - (e * sin(L_b))^2);

% Convert position using (2.112)

cos_lat = cos(L_b);
sin_lat = sin(L_b);
cos_long = cos(lambda_b);
sin_long = sin(lambda_b);
x = (R_E + h_b) * cos_lat * cos_long;
y = (R_E + h_b) * cos_lat * sin_long;
z = ((1 - e^2) * R_E + h_b) * sin_lat;
pos = [x,y,z];
          
          