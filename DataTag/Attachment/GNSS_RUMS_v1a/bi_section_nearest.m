function M_xyz = bi_section_nearest(c1_xyz,c2_xyz,rcvr_xyz,sv_xyz)

% find the nearest point along line of c1-c2 between SV & RCVR
% para_reso control the resolution on edge finding the point

para_reso = 0.019/10;
% D2R = pi/180;
% R2D = 180/pi;

A_xyz = c1_xyz;
B_xyz = c2_xyz;
% A_xyz = llh2xyz([22.299544,114.178556,52].*[D2R,D2R,1]);
% A_xyz0 = A_xyz;
M_xyz = (A_xyz+B_xyz)./2;
sec_reso = norm(A_xyz - B_xyz)/2;

while sec_reso > para_reso
    A_dist = norm(sv_xyz - A_xyz) + norm(A_xyz - rcvr_xyz) - norm(sv_xyz - rcvr_xyz);
    B_dist = norm(sv_xyz - B_xyz) + norm(B_xyz - rcvr_xyz) - norm(sv_xyz - rcvr_xyz);
%     M_dist = norm(sv_xyz - M_xyz) + norm(M_xyz - rcvr_xyz) - norm(sv_xyz - rcvr_xyz);
    if B_dist > A_dist
        B_xyz = M_xyz;
    else
        A_xyz = M_xyz;
    end
    M_xyz = (A_xyz+B_xyz)./2;
    sec_reso = norm(A_xyz - B_xyz)/2;
end 
    
% M_llh = xyz2llh(M_xyz);
% M_llh = M_llh.*[R2D,R2D,1];
%     
% 
% norm(A_xyz0-M_xyz)

% A_xyz = c1_xyz;
% B_xyz = c2_xyz;
% A_xyz = llh2xyz([22.299544,114.178556,52].*[D2R,D2R,1]);
%     for ratio = 0:0.01:1
%         eid = round(ratio*100+1);
%         inter_xyz(eid,1:3) = A_xyz+(B_xyz-A_xyz).*ratio;
%         [inter_elaz(eid,2),inter_elaz(eid,1),~] = topocent(rcvr_xyz,inter_xyz(eid,1:3));
%         inter_dist(eid,1) = norm(inter_xyz(eid,1:3)-sv_xyz)+norm(inter_xyz(eid,1:3)-rcvr_xyz) - norm(sv_xyz-rcvr_xyz);
%     end
%     [~,eid] = min(inter_dist);
%     diff_xyz = A_xyz+(B_xyz-A_xyz).*(eid-1)/100;
%     M_llh = xyz2llh(diff_xyz);
% M_llh = M_llh.*[R2D,R2D,1];



    