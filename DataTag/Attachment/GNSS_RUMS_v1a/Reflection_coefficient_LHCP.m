function [R_LHCP] = Reflection_coefficient_LHCP(reflecting_pos,SV_pos,user_pos,lambda)

% calculating LHCP signal reflection coefficient
% Referenced by Taro suggested paper
%
% Jia, Yan, and Yuekun Pei. "Remote Sensing in Land Applications by Using 
% GNSS-Reflectometry." Recent Advances and Applications in Remote Sensing. 
% IntechOpen, 2018.

R_LHCP = zeros(size(reflecting_pos,1),1);
for idn = 1:1:size(reflecting_pos,1)
        c = norm(SV_pos-user_pos);
        a = norm(SV_pos-reflecting_pos(idn,:));
        b = norm(user_pos-reflecting_pos(idn,:));
        i_angle = acosd((a^2+b^2-c^2)/(2*a*b))/2;%incidence angle
%         ci_angle = 90-i_angle;%complementary angle of incidence
        
        % ---Taro method ---
        % refractive index 
%         n1 = 1.000293;%Air
%         n2 = 1.52;%Window glass
%         n12 = n2/n1;
        % magnetic permeability
%         miu1 = 1.25663753*10^-6;%Air
%         miu2 = 4*pi*10^-7;%Concrete,glass
        % perpendicular polarized electric field
%         R_v = (miu2*sind(ci_angle)-miu1*sqrt(n12^2-cosd(ci_angle)^2))/...
%               (miu2*sind(ci_angle)+miu1*sqrt(n12^2-cosd(ci_angle)^2));
        % parallel polarized electric field
%         R_h = (miu1*n12^2*sind(ci_angle)-miu2*sqrt(n12^2-cosd(ci_angle)^2))/...
%               (miu1*n12^2*sind(ci_angle)+miu2*sqrt(n12^2-cosd(ci_angle)^2));
        % LHCP reflection coefficient
%         R_LHCP_taro = (R_v-R_h)/2;
        
        % ---Yan method---
        % dielectric constant
        ep = 4.7;%Glass(4.7)/Concrete(4.5)
        ep0 = 1.00058986;%Air
        % electrical conductivity
%         cigma = 10^-11;
%         rho = 10^-8;
%         cigma = 1/rho;
%         ep_r = ep/ep0-60j*lambda*cigma;
        ep_r = ep/ep0; % lossy medium conduction term neglected
        % horizontal polarized electric field
        R_h = (cosd(i_angle)-sqrt(ep_r-sind(i_angle)^2))/...
              (cosd(i_angle)+sqrt(ep_r-sind(i_angle)^2));
        % vertical polarized electric field
        R_v = (ep_r*cosd(i_angle)-sqrt(ep_r-sind(i_angle)^2))/...
              (ep_r*cosd(i_angle)+sqrt(ep_r-sind(i_angle)^2));
        % LHCP reflection coefficient
        R_LHCP(idn,1) = (R_v-R_h)/2;
%         disp(['R_LHCP: ',num2str(R_LHCP_taro),'(Taro); ',num2str(R_LHCP_yan),'(Yan)']);
end

