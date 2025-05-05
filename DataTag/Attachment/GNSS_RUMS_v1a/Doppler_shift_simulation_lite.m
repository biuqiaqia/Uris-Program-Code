function [dop_LOS,dop_R] = Doppler_shift_simulation_lite(XR,VR,XS,VS,dtRV,lambda,vis_tag,inter_pos)
% ------------------------------------------------------------------------
% Simulate GNSS Doppler shift measurements
%
% by GH.Zhang 2020/09/08
% guo-hao.zhang@connect.polyu.hk
% ------------------------------------------------------------------------
% Input:
%       XR - User position in ECEF x,y,z (m)
%       VR - User velocity in ECEF x,y,z (m/s)
%       XS - Satellite positions in ECEF x,y,z (m)
%       VS - Satellite velocities in ECEF x,y,z (m)
%     dtRV - Receiver clock drift (m/s)
%   lambda - Wavelength of each signal (m)
%  vis_tag - Visibility tag (1-LOS,2-MP,3-Reflection,4-Diffraction,5-R&D)
%
% Output:
%    dop_R - simulated Doppler shift (Hz)
% ------------------------------------------------------------------------

v_light = goGNSS.V_LIGHT;
dtRV = -dtRV;%LOS-vector oppose between GoGPS/Kaplan

% LOS Doppler
a = (XS' - XR')./norm(XS - XR);%User-SV unit vector
dop_LOS = (VR*a-VS*a-v_light*dtRV)/lambda;

% type-based Doppler
switch vis_tag
    case 0
        dop_R = nan;
    case 1 %LOS
        a = (XS' - XR')./norm(XS - XR);%User-SV unit vector
        dop_R = (VR*a-VS*a-v_light*dtRV)/lambda;
	case 2 %intersection
        a1 = (XS' - inter_pos')./norm(XS - inter_pos);%Interpoint-SV unit vector
        a2 = (inter_pos' - XR')./norm(inter_pos - XR);%User-Interpoint unit vector
        dop_R = (VR*a2-VS*a1-v_light*dtRV-(VR*a2)*(VS*a1)/v_light)/lambda;
end







% switch vis_tag
%     case 0
%         dop_R = nan;
%     case 1 %LOS
%         a = (XS' - XR')./norm(XS - XR);%User-SV unit vector
%         dop_R = (VR*a-VS*a-v_light*dtRV)/lambda;
% 	case 2 %MP
%         a = (XS' - XR')./norm(XS - XR);%User-SV unit vector
%         dop_R = (VR*a-VS*a-v_light*dtRV)/lambda;
% 	case 3 %Reflection
%         a1 = (XS' - inter_pos')./norm(XS - inter_pos);%Interpoint-SV unit vector
%         a2 = (inter_pos' - XR')./norm(inter_pos - XR);%User-Interpoint unit vector
%         dop_R = (VR*a2-VS*a1-v_light*dtRV-(VR*a2)*(VS*a1)/v_light)/lambda;
%     case 4 %Diffraction
%         a1 = (XS' - inter_pos')./norm(XS - inter_pos);%Interpoint-SV unit vector
%         a2 = (inter_pos' - XR')./norm(inter_pos - XR);%User-Interpoint unit vector
%         dop_R = (VR*a2-VS*a1-v_light*dtRV-(VR*a2)*(VS*a1)/v_light)/lambda;
%     case 5 %R&D
%         a1 = (XS' - inter_pos')./norm(XS - inter_pos);%Interpoint-SV unit vector
%         a2 = (inter_pos' - XR')./norm(inter_pos - XR);%User-Interpoint unit vector
%         dop_R = (VR*a2-VS*a1-v_light*dtRV-(VR*a2)*(VS*a1)/v_light)/lambda;
% %         dop_R = nan;
% end





%number of observations
% n = size(XS,1);
% %number of unknown parameters
% m = 4;
% %approximate receiver-satellite distance
% XR_mat = XR(ones(1,1),:);
% distR_approx = sqrt(sum((XS-XR_mat).^2 ,2));
% %geometry matrix (same as Kaplan's book)
% A = [(XR(1) - XS(:,1)) ./ distR_approx, ... %column for X coordinate
%      (XR(2) - XS(:,2)) ./ distR_approx, ... %column for Y coordinate
%      (XR(3) - XS(:,3)) ./ distR_approx, ... %column for Z coordinate
%      ones(1,1)];        %column for receiver clock delay (multiplied by c)
% 
% % QZSS clk equals to GPS
% sys(find(sys==5)) = 1;
% uni_sys = unique(sys(sys ~= 0));
% num_sys = length(uni_sys);
% ISB = zeros(n,1);
% if (num_sys > 1)
%     m = m + num_sys - 1;
%     for s = 2 : num_sys
%         ISB(sys == uni_sys(s)) = 1;
%         A = [A, ISB];
%         ISB = zeros(n,1);
%     end
% end
% %known term vector
% b = sum(A(:,1:3).*VS,2);
% x = [-VR,dtRV*v_light]';
% y0 = A*x + b;
% dop_R = y0./lambda';