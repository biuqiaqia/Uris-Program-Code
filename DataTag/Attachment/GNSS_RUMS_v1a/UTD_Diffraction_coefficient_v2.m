function [D_total] = UTD_Diffraction_coefficient_v2(beta0,phi,phi_,s,lambda)

% Calculating diffraction coefficient based on UTD model
% Version 2 verified by knife-edge model
% 
% Reference:
% [1] Taro Suzuki and Nobuaki Kubo. "GNSS positioning with multipath 
%     simulation using 3D surface model in urban canyon." Proceedings of 
%     the 25th International Technical Meeting of the Satellite Division 
%     of The Institute of Navigation (ION GNSS 2012). 2001.
% [2] Jia, Yan, and Yuekun Pei. "Remote Sensing in Land Applications by Using 
%     GNSS-Reflectometry." Recent Advances and Applications in Remote Sensing. 
%     IntechOpen, 2018.
% [3] McNamara, D. A., C. W. I. Pistorius, and J. A. G. Malherbe. 
%     "Introduction to The Uniform Geometrical theory of diffraction." 
%     Artech House, London (1990).
% [4] Kouyoumjian, R. G. and P. H. Pathak (1974). "A uniform geometrical 
%     theory of diffraction for an edge in a perfectly conducting surface." 
%     Proceedings of the IEEE 62(11): 1448-1461.


% constant parameters
D2R = pi/180;
% R2D = 180/pi;
% interior wedge parameter based on interior wedge angle alpha (rad)
alpha = pi/2;%building edge
n = (2*pi-alpha)/pi;
% wave number
k = 2*pi/lambda;

% cone half-angle
% beta0 = beta0;
% diffracted-edge-perpendicular angle
% phi = phi;
% incident-edge-perpendicular angle
% phi_ = phi_;
% distance parameter with diffracted point 
L = s*sind(beta0)*sind(beta0);%(plan-wave approximation)

% a+ & a- function
miu_n = (phi - phi_)*D2R;
% miu_p = (phi + phi_)*D2R;
N_p_id = 0;
for N_p = -2:1:2
    N_p_id = N_p_id + 1;
    delta_p(N_p_id,1) = N_p;
    delta_p(N_p_id,2) = abs(2*pi*n*N_p - miu_n - pi);
end
[~,N_p_id] = min(delta_p(:,2));
N_p = delta_p(N_p_id,1);%determine N+
N_n_id = 0;
for N_n = -2:1:2
    N_n_id = N_n_id + 1;
    delta_n(N_n_id,1) = N_n;
    delta_n(N_n_id,2) = abs(2*pi*n*N_n - miu_n + pi);
end
[~,N_n_id] = min(delta_n(:,2));
N_n = delta_n(N_n_id,1);%determine N-
a_p_n = 2*cos((2*n*pi*N_p-miu_n)/2)^2;
a_n_n = 2*cos((2*n*pi*N_n-miu_n)/2)^2;
% a_p_p = 2*cos((2*n*pi*N_p-miu_p)/2)^2;
% a_n_p = 2*cos((2*n*pi*N_n-miu_p)/2)^2;

% Fresnel integral
F_p_n = Fresnel_integral(k*L*a_p_n);
F_n_n = Fresnel_integral(k*L*a_n_n);
% F_p_p = Fresnel_integral(k*Lrn*a_p_p);
% F_n_p = Fresnel_integral(k*Lro*a_n_p);

% UTD model
D1 = (-exp(-1i*pi/4)/(2*n*sqrt(2*pi*k)*sind(beta0)))*(cot((pi+miu_n)/2/n)*F_p_n);
D2 = (-exp(-1i*pi/4)/(2*n*sqrt(2*pi*k)*sind(beta0)))*(cot((pi-miu_n)/2/n)*F_n_n);
% D3 = (-exp(-1i*pi/4)/(2*n*sqrt(2*pi*k)*sind(beta0)))*(cot((pi+miu_p)/2/n)*F_p_p);
% D4 = (-exp(-1i*pi/4)/(2*n*sqrt(2*pi*k)*sind(beta0)))*(cot((pi-miu_p)/2/n)*F_n_p);
D_total = [D1+D2,D1,D2];
















% for each possible diffraction
% for idn = 1:1:size(diffraction_set{1,1},1)
% diff_xyz = diffraction_set{1,1}(idn,:);
% diff_corner = diffraction_set{1,3}(idn,1:3);
% diff_enu = xyz2enu(diff_xyz,user_pos);
% diff_c_enu = xyz2enu(diff_corner,user_pos);
% sv_enu = xyz2enu(SV_pos,user_pos);

% incident-edge angle
% a = norm(diff_xyz-SV_pos);
% b = norm(SV_pos-diff_corner);
% c = norm(diff_xyz-diff_corner);
% beta_ = acosd((a^2+c^2-b^2)/(2*a*c));
% i_tag = 1;
% if beta_>90
%     beta_ = 180 - beta_;
%     i_tag = -1;
% end

% diffracted-edge angle
% a = norm(diff_xyz-user_pos);
% b = norm(user_pos-diff_corner);
% c = norm(diff_xyz-diff_corner);
% beta = acosd((a^2+c^2-b^2)/(2*a*c));
% d_tag = 1;
% if beta>90
%     beta = 180 - beta;
%     d_tag = -1;
% end

% diffracted-edge-perpendicular angle [ENU]
% v_base = diff_c_enu - diff_enu;
% v_diff = - diff_enu;
% d_diff2base = norm(v_diff'*v_base)/norm(v_base);
% r_diff2base = d_diff2base/norm(v_base)*d_tag;
% p_diff_enu = diff_enu + r_diff2base * v_base;
% inter_angel = cosd(p_diff_enu(3)/norm(p_diff_enu));
% phi = 270-inter_angel;

% incident-edge-perpendicular angle [ENU]
% v_inci = sv_enu - diff_enu;
% d_inci2base = norm(v_inci'*v_base)/norm(v_base);
% r_inci2base = d_inci2base/norm(v_base)*i_tag;
% p_inci_enu = diff_enu + r_inci2base * v_base;
% dp_in = norm(sv_enu-p_inci_enu);
% dp_h = p_inci_enu(3);
% dp_l = norm(sv_enu-[p_inci_enu(1);p_inci_enu(2);0]);
% phi_ = acosd((dp_in^2+dp_h^2-dp_l^2)/(2*dp_in*dp_h)) + inter_angel;

% output angle [degree]
% diff_angle(idn,1:4) = [beta_,beta,phi_,phi];

% distance parameter
% s = norm(user_pos-diff_xyz);
% s_= norm(SV_pos-diff_xyz);
% L = (s*s_)/(s+s_)*sind(beta)^2;

% a+ & a- function
% miu_n = (phi - phi_)*D2R;
% miu_p = (phi + phi_)*D2R;
% N_p_id = 0;
% for N_p = -2:1:2
%     N_p_id = N_p_id + 1;
%     delta_p(N_p_id,1) = N_p;
%     delta_p(N_p_id,2) = abs(2*pi*n*N_p - miu_n - pi);
% end
% [~,N_p_id] = min(delta_p(:,2));
% N_p = delta_p(N_p_id,1);%determine N+
% N_n_id = 0;
% for N_n = -2:1:2
%     N_n_id = N_n_id + 1;
%     delta_n(N_n_id,1) = N_n;
%     delta_n(N_n_id,2) = abs(2*pi*n*N_n - miu_n + pi);
% end
% [~,N_n_id] = min(delta_n(:,2));
% N_n = delta_n(N_n_id,1);%determine N-
% a_p_n = 2*cos((2*n*pi*N_p-miu_n)/2)^2;
% a_n_n = 2*cos((2*n*pi*N_n-miu_n)/2)^2;
% a_p_p = 2*cos((2*n*pi*N_p-miu_p)/2)^2;
% a_n_p = 2*cos((2*n*pi*N_n-miu_p)/2)^2;
% 
% % Fresnel integral
% F_p_n = Fresnel_integral(k*L*a_p_n);
% F_n_n = Fresnel_integral(k*L*a_n_n);
% F_p_p = Fresnel_integral(k*L*a_p_p);
% F_n_p = Fresnel_integral(k*L*a_n_p);
% 
% % UTD model
% % [~,idd] = min(diffraction_set{1,2});
% % R_reflect = R_reflect_set(idd,1);
% % R_reflect = 0;
% D1 = (-exp(-1i*pi/4)/(2*n*sqrt(2*pi*k)*sind(beta0)))*(cot((pi+miu_n)/2/n)*F_p_n);
% D2 = (-exp(-1i*pi/4)/(2*n*sqrt(2*pi*k)*sind(beta0)))*(cot((pi-miu_n)/2/n)*F_n_n);
% D3 = (-exp(-1i*pi/4)/(2*n*sqrt(2*pi*k)*sind(beta0)))*(cot((pi+miu_p)/2/n)*F_p_p);
% D4 = (-exp(-1i*pi/4)/(2*n*sqrt(2*pi*k)*sind(beta0)))*(cot((pi-miu_p)/2/n)*F_n_p);


% D_i = (-exp(-1i*pi/4)/(2*n*sqrt(2*pi*k)*sind(beta)))*...
%       (cot((pi+miu_n)/2/n)*F_p_n + R_reflect*cot((pi-miu_p)/2/n)*F_n_n);
% D_c = (-exp(-1i*pi/4)/(2*n*sqrt(2*pi*k)*sind(beta)))*...
%       (cot((pi+miu_n)/2/n)*F_p_n + R_reflect*cot((pi-miu_p)/2/n)*F_n_p);
% D_all = (-exp(-1i*pi/4)/(2*n*sqrt(2*pi*k)*sind(beta)))*...
%         (cot((pi+miu_n)/2/n)*F_p_n + ...
%          cot((pi-miu_n)/2/n)*F_n_n + ...
%          R_reflect*cot((pi+miu_p)/2/n)*F_p_p + ...
%          R_reflect*cot((pi-miu_p)/2/n)*F_n_p);
  
% D_coeff(idn,1:3) = [D_i,D_c,D_all];
% D_coeff(idn,1) = D_all;

% Individal electric field
% E_i_ratio(idn,1) = sqrt(1/(s*s_*(s+s_)))*D_i*exp(-1i*k*(s+s_));% from taro

% E_i_ratio(idn,1) = sqrt(s_/(s*(s+s_)))*D_all*exp(-1i*k*s);

% check_E_ratio(idn,1) = sqrt(s_/(s*(s+s_)))*D_i*exp(-1i*k*s);% revised from [4]
% check_E_ratio(idn,2) = sqrt(s_/(s*(s+s_)))*D_c*exp(-1i*k*s);
% check_E_ratio(idn,3) = sqrt(s_/(s*(s+s_)))*D_all*exp(-1i*k*s);

% end














% % incident angle related to edge
% beta = nan*D2R;%---need update
% % incident angle related to edge-perpendicular plane
% phi = nan*D2R;
% % diffracted angle related to edge-perpendicular plane
% phi_ = nan*D2R;
% 

% 
% % 
% 
% 

% 
% 
% 
% 
% 
% 
% 
% 
% R_LHCP = zeros(size(reflecting_pos,1),1);
% for idn = 1:1:size(reflecting_pos,1)
%         c = norm(SV_pos-user_pos);
%         a = norm(SV_pos-reflecting_pos(idn,:));
%         b = norm(user_pos-reflecting_pos(idn,:));
%         i_angle = acosd((a^2+b^2-c^2)/(2*a*b))/2;%incidence angle
% %         ci_angle = 90-i_angle;%complementary angle of incidence
%         
%         % ---Taro method ---
%         % refractive index 
% %         n1 = 1.000293;%Air
% %         n2 = 1.52;%Window glass
% %         n12 = n2/n1;
%         % magnetic permeability
% %         miu1 = 1.25663753*10^-6;%Air
% %         miu2 = 4*pi*10^-7;%Concrete,glass
%         % perpendicular polarized electric field
% %         R_v = (miu2*sind(ci_angle)-miu1*sqrt(n12^2-cosd(ci_angle)^2))/...
% %               (miu2*sind(ci_angle)+miu1*sqrt(n12^2-cosd(ci_angle)^2));
%         % parallel polarized electric field
% %         R_h = (miu1*n12^2*sind(ci_angle)-miu2*sqrt(n12^2-cosd(ci_angle)^2))/...
% %               (miu1*n12^2*sind(ci_angle)+miu2*sqrt(n12^2-cosd(ci_angle)^2));
%         % LHCP reflection coefficient
% %         R_LHCP_taro = (R_v-R_h)/2;
%         
%         % ---Yan method---
%         % dielectric constant
%         ep = 4.7;%Glass(4.7)/Concrete(4.5)
%         ep0 = 1.00058986;%Air
%         % electrical conductivity
% %         cigma = 10^-11;
% %         ep_r = ep/ep0-60j*lambda(ids)*cigma*0;
%         ep_r = ep/ep0;
%         % horizontal polarized electric field
%         R_h = (cosd(i_angle)-sqrt(ep_r-sind(i_angle)^2))/...
%               (cosd(i_angle)+sqrt(ep_r-sind(i_angle)^2));
%         % vertical polarized electric field
%         R_v = (ep_r*cosd(i_angle)-sqrt(ep_r-sind(i_angle)^2))/...
%               (ep_r*cosd(i_angle)+sqrt(ep_r-sind(i_angle)^2));
%         % LHCP reflection coefficient
%         R_LHCP(idn,1) = (R_v-R_h)/2;
% %         disp(['R_LHCP: ',num2str(R_LHCP_taro),'(Taro); ',num2str(R_LHCP_yan),'(Yan)']);
% end

