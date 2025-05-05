function [pos_est,pr_resi,DOP,dtR_est] = Least_square_positioning(sv_data,pos_ini)

% -------------------------------------------------------------------------
% Least square positioning
% Multi-constellation (GPS/GLONASS/GALILEO/BEIDOU)
% Pseudorange residual and DOP supported
% Input:    sv_data =   column 1: PRN code
%                       column 2: Pseudorange
%                       column 3-5: Satellite position 
%           pos_ini =   initial guess of position in ECEF
%
% Output:   pos_est =   position estimation
%           pr_resi =   pseudorange residual
%           DOP     =   HDOP,VDOP
%
% by GH.Zhang 2019/03/07
% -------------------------------------------------------------------------
H_t=[];
constellation = zeros(size(sv_data,1),1);
for idm = 1:1:size(sv_data,1)
    if sv_data(idm,1) <= 32
        constellation(idm) = 1; 
        H_t(idm,:) = [1,0,0,0];
    elseif sv_data(idm,1) > 32 && sv_data(idm,1) <= 56
        constellation(idm) = 2;
        H_t(idm,:) = [1,1,0,0];
    elseif sv_data(idm,1) > 56 && sv_data(idm,1) <= 86
        constellation(idm) = 3;
        H_t(idm,:) = [1,0,1,0];
    elseif sv_data(idm,1) > 86 && sv_data(idm,1) <= 123
        constellation(idm) = 4;
        H_t(idm,:) = [1,0,0,1];
    elseif sv_data(idm,1) > 123 && sv_data(idm,1) <= 127
        constellation(idm) = 1;
        H_t(idm,:) = [1,0,0,0];
    end
end
H_t = H_t(:,ismember([1;2;3;4],constellation));

if size(sv_data,1)>=(size(H_t,2)+3)
    dt = zeros(size(H_t,2),1);
    delta_pos = [1e9; 1e9; 1e9; dt];

    while norm(delta_pos(1:3,1))>1e-3
        clear H;
        for id_sv = 1:1:size(sv_data,1)
            pos_pr0(id_sv) = norm(sv_data(id_sv,3:5)-pos_ini);
            y(id_sv,1) = sv_data(id_sv,2) - pos_pr0(id_sv) - dt(1,1);
            H_v(id_sv,1:3) = (pos_ini - sv_data(id_sv,3:5))./pos_pr0(id_sv);
            %DOP
            sv_enu = xyz2enu(sv_data(id_sv,3:5),pos_ini);
            H_dop(id_sv,1:3) = sv_enu'./pos_pr0(id_sv);
        end
        H = [H_v,H_t];
        delta_pos = inv(H'*H)*H'*y;
        pos_ini = pos_ini + delta_pos(1:3,1)';
        dt = dt*1 + delta_pos(4:end,1);
    end
    pos_est = pos_ini;
    pr_resi = y - H*delta_pos;%residual
    H_dop = [H_dop,H_t];
    Q = (H_dop'*H_dop)^-1;
    DOP = [sqrt(Q(1,1)+Q(2,2)),sqrt(Q(3,3))];%HDOP,VDOP
    dtR_est = H_t*dt;
%     dtR_est = dt;
else
    pos_est = [nan,nan,nan];
    pr_resi = nan(size(sv_data,1),1);
    DOP = [nan,nan];
    dtR_est = nan;
end



















