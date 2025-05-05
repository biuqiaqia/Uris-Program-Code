function [HDOP,VDOP] = DOP_cal(sv_pos,pos_ini)
% calculate the DOP from satellite position in ECEF

% sv_pos = data{2,2}{1,3}(:,3:5);
% pos_ini = data{2,2}{1,2};

for idsv = 1:1:size(sv_pos,1)
    R(idsv,1) = norm(sv_pos(idsv,1:3)-pos_ini);
    sv_enu = xyz2enu(sv_pos(idsv,1:3),pos_ini);
    H(idsv,1:4) = [sv_enu'./R(idsv,1),1];
end

Q = (H'*H)^-1;
HDOP = sqrt(Q(1,1)+Q(2,2));
VDOP = sqrt(Q(3,3));



