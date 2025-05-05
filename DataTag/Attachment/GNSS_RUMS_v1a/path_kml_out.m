function inc_ang = path_kml_out(inter_pos,sv_pos,rcvr_pos,color,out_path)

% D2R = pi/180;
R2D = 180/pi;

sv_llh = xyz2llh(sv_pos);
sv_llh = sv_llh.*[R2D,R2D,1];
rcvr_llh = xyz2llh(rcvr_pos);
rcvr_llh = rcvr_llh.*[R2D,R2D,1];
inter_llh = xyz2llh(inter_pos);
inter_llh = inter_llh.*[R2D,R2D,1];

if isempty(color)
    color = [1,0,0];
end

s = norm(inter_pos-rcvr_pos);
s_ = norm(inter_pos-sv_pos);
s0 = norm(rcvr_pos-sv_pos);
inc_ang = acosd((s^2+s_^2-s0^2)/(2*s*s_))/2;

kmlwriteline(out_path,...
                 [sv_llh(1),inter_llh(1),rcvr_llh(1)],...
                 [sv_llh(2),inter_llh(2),rcvr_llh(2)],...
                 [sv_llh(3),inter_llh(3),rcvr_llh(3)],...
                 'Color',color,'LineWidth',5,'Name','  ');

