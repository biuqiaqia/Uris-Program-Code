function [cnr_diff,D_UTD] = UTD_GNSS_CNR_sim(rcvr_xyz,sv_xyz,diff_xyz,corner_pos,inter_coeff,inter_delay,lambda,cnr_open)

% ------------------------------------------------------------------------
% The Uniform Theory of Diffraction simulation on GNSS
% Diffracted Carrier-to-noise ratio
%
% by GH.Zhang 2020/07/14
% guo-hao.zhang@connect.polyu.hk
% ------------------------------------------------------------------------

% [sv_az,sv_el,sv_dist] = topocent(rcvr_xyz,sv_xyz);
s = norm(diff_xyz-rcvr_xyz);
s_ = norm(diff_xyz-sv_xyz);
s0 = norm(rcvr_xyz-sv_xyz);
diff_delay = s_ + s - s0;

% o-face corner_xyz
c1_xyz = corner_pos(1,1:3);
c2_xyz = corner_pos(2,1:3);
c3_xyz = corner_pos(3,1:3);

v_s_prime = (diff_xyz-sv_xyz)./norm(diff_xyz-sv_xyz);
v_s = (rcvr_xyz-diff_xyz)./norm(rcvr_xyz-diff_xyz);
v_e = (c2_xyz-c1_xyz)./norm(c2_xyz-c1_xyz);
v_3 = (c3_xyz-c1_xyz)./norm(c3_xyz-c1_xyz);
v_n_prime = cross(v_e,v_3)./norm(cross(v_e,v_3));%new one

v_phi_prime = -cross(v_e,v_s_prime)./norm(cross(v_e,v_s_prime));
v_phi = cross(v_e,v_s)./norm(cross(v_e,v_s));
v_beta0_prime = cross(v_phi_prime,v_s_prime);
v_beta0 = cross(v_phi,v_s);

beta0 = acosd(dot(v_e,v_s_prime));
phi_ = acosd(dot(v_n_prime,v_phi_prime));
phi = 360-acosd(dot(v_n_prime,v_phi));%angle overflow

[D3_sh(1,1:3)] = UTD_Diffraction_coefficient_v2(beta0,phi,phi_,s,lambda);% [D1+D2,D1,D2]
D3_RR = -dot(v_beta0_prime,v_beta0)*D3_sh(1,1)*0.5-dot(v_phi_prime,v_phi)*D3_sh(1,1)*0.5;% dyadic form for RHCP-RHCP

% interference between LOS/Diff
% if vis_tag == 0
%     D_UTD_diff = D3_RR*sqrt(1/s)*exp(-1i*(2*pi/lambda)*s);
% 	D_UTD_all = D3_RR*sqrt(1/s)*exp(-1i*(2*pi/lambda)*s);
% elseif vis_tag == 1
%     D_UTD_diff = D3_RR*sqrt(1/s)*exp(-1i*(2*pi/lambda)*s);
% 	D_UTD_all = D3_RR*sqrt(1/s)*exp(-1i*(2*pi/lambda)*s) + exp(-1i*(2*pi/lambda)*(s-diff_delay));
% end

% interference between Diff with other field with coefficient
% (0-no interference/1-LOS/R-reflected field)
D_UTD_diff = D3_RR*sqrt(1/s)*exp(-1i*(2*pi/lambda)*s);
if length(inter_coeff)==1
    D_UTD_all = D3_RR*sqrt(1/s)*exp(-1i*(2*pi/lambda)*s) + inter_coeff*exp(-1i*(2*pi/lambda)*(s-diff_delay+inter_delay));
elseif length(inter_coeff)==2
    D_UTD_all = D3_RR*sqrt(1/s)*exp(-1i*(2*pi/lambda)*s) + ...
                inter_coeff(1)*exp(-1i*(2*pi/lambda)*(s-diff_delay+inter_delay(1))) + ...
                inter_coeff(2)*exp(-1i*(2*pi/lambda)*(s-diff_delay+inter_delay(2)));
end
D_UTD = [D_UTD_all,D_UTD_diff,D3_RR,D3_sh];

% CNR attenuation by diffraction coefficient
cnr_diff = cnr_open + 10*log10((s0/(s_+s)*abs(D_UTD_all))^2);

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
