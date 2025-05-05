function [pos_est,pr_resi,tag] = LS_mpos_pr_resi(pr,sv_pos,pos_ini)

% Least square positioning (conventional with multi-constellation)
% Input:    pr      =   pseudorange measurement
%           sv_pos  =   satellite position in ECEF
%           pos_ini =   initial guess of position in ECEF
%
% Output:   pos_est =   position estimation
%           dt      =   receiver clock bias or SSE
%
% modified in 2018/03/29 by GH.ZHANG

% pr = pr_check;
% sv_pos = SV;
% pos_ini = c_data(1,3:5);

%% Constellation tag
tag_N = 0;
tag_G = 0;
tag_B = 0;
for idc = 1:1:size(sv_pos,1)
    if sv_pos(idc,1) <= 32
        tag_N = 1;
    elseif sv_pos(idc,1) >= 56
        tag_B = 1;
    elseif (sv_pos(idc,1) > 32) && (sv_pos(idc,1) <= 56)
        tag_G = 1;
    end
end

%% Absolute positioning with conventional LS
if size(pr,1)<4
    pos_est = 0;
elseif tag_N && ~tag_G && ~tag_B %GPS only
    [pos_est,pr_resi] = LS_pos_N(pr,sv_pos,pos_ini);
    tag = 1;
elseif ~tag_N && tag_G && ~tag_B %GLONASS only
    [pos_est,pr_resi] = LS_pos_G(pr,sv_pos,pos_ini);
    tag = 2;
elseif ~tag_N && ~tag_G && tag_B %BDS only
    [pos_est,pr_resi] = LS_pos_B(pr,sv_pos,pos_ini);
    tag = 3;
elseif tag_N && tag_G && ~tag_B %GPS/GLONASS
    [pos_est,pr_resi] = LS_pos_NG(pr,sv_pos,pos_ini);
    tag = 4;
elseif tag_N && ~tag_G && tag_B %GPS/BDS
    [pos_est,pr_resi] = LS_pos_NB(pr,sv_pos,pos_ini);
    tag = 5;
 elseif ~tag_N && tag_G && tag_B %BDS/GLONASS
    [pos_est,pr_resi] = LS_pos_GB(pr,sv_pos,pos_ini);
    tag = 6;
elseif tag_N && tag_G && tag_B %GPS/GLONASS/BDS
    [pos_est,pr_resi] = LS_pos_NGB(pr,sv_pos,pos_ini);
    tag = 7;
end
function [pos_est,pr_resi] = LS_pos_N(pr,sv_pos,pos_ini)

% Least square positioning (GPS only)
% Input:    pr      =   pseudorange measurement
%           sv_pos  =   satellite position in ECEF
%           pos_ini =   initial guess of position in ECEF
%
% Output:   pos_est =   position estimation
%           sse     =   sum of square error

if size(pr,1)>=4
    dt = [0,0,0];
    delta_pos = [1e9; 1e9; 1e9; 0];
    while ((norm(delta_pos(1:3,1))>1e-3))
        for idx_sv = 1 : size(sv_pos,1)
            pos1_pr0(idx_sv) = norm(sv_pos(idx_sv,2:4)-pos_ini);
            prn = sv_pos(idx_sv,1);
            if prn <= 32
                y(idx_sv,1) = pr(idx_sv) - pos1_pr0(idx_sv) - delta_pos(4,1);
                H(idx_sv,1:4) = [(pos_ini - sv_pos(idx_sv,2:4))./pos1_pr0(idx_sv),1];
            end
        end
    delta_pos = inv(H'*H)*H'*y;
    pos_ini = pos_ini + delta_pos(1:3,1)';
    dt(1,1) = dt(1,1) + delta_pos(4,1);
    end
    pos_est = pos_ini;
    pr_resi = y - H*delta_pos;%residual
    tag = 1;
end

function [pos_est,pr_resi] = LS_pos_B(pr,sv_pos,pos_ini)
if size(pr,1)>=4
    dt = [0,0,0];
    delta_pos = [1e9; 1e9; 1e9; 0];
    while ((norm(delta_pos(1:3,1))>1e-3))
        for idx_sv = 1 : size(sv_pos,1)
            pos1_pr0(idx_sv) = norm(sv_pos(idx_sv,2:4)-pos_ini);
            prn = sv_pos(idx_sv,1);
            if prn > 86
                y(idx_sv,1) = pr(idx_sv) - pos1_pr0(idx_sv) - delta_pos(4,1);
                H(idx_sv,1:4) = [(pos_ini - sv_pos(idx_sv,2:4))./pos1_pr0(idx_sv),1];
            end
        end
    delta_pos = inv(H'*H)*H'*y;
    pos_ini = pos_ini + delta_pos(1:3,1)';
    dt(1,3) = dt(1,3) + delta_pos(4,1);
    end
    pos_est = pos_ini;
    pr_resi = y - H*delta_pos;
    tag = 3;
end

function [pos_est,pr_resi] = LS_pos_G(pr,sv_pos,pos_ini)

% Least square positioning (GLONASS only)
% Input:    pr      =   pseudorange measurement
%           sv_pos  =   satellite position in ECEF
%           pos_ini =   initial guess of position in ECEF
%
% Output:   pos_est =   position estimation
%           sse     =   sum of square error

if size(pr,1)>=4
    dt = [0,0,0];
    delta_pos = [1e9; 1e9; 1e9; 0];
    while ((norm(delta_pos(1:3,1))>1e-3))
        for idx_sv = 1 : size(sv_pos,1)
            pos1_pr0(idx_sv) = norm(sv_pos(idx_sv,2:4)-pos_ini);
            prn = sv_pos(idx_sv,1);
            if prn > 32 && prn <= 56
                y(idx_sv,1) = pr(idx_sv) - pos1_pr0(idx_sv) - delta_pos(4,1);
                H(idx_sv,1:4) = [(pos_ini - sv_pos(idx_sv,2:4))./pos1_pr0(idx_sv),1];
            end
        end
    delta_pos = inv(H'*H)*H'*y;
    pos_ini = pos_ini + delta_pos(1:3,1)';
    dt(1,2) = dt(1,2) + delta_pos(4,1);
    end
    pos_est = pos_ini;
    pr_resi = y - H*delta_pos;
    tag = 2;
end

function [pos_est,pr_resi] = LS_pos_GB(pr,sv_pos,pos_ini)

% Least square positioning (GLONASS & BDS)
% Input:    pr      =   pseudorange measurement
%           sv_pos  =   satellite position in ECEF
%           pos_ini =   initial guess of position in ECEF
%
% Output:   pos_est =   position estimation
%           sse     =   sum of square error

if size(pr,1)>=5
    dt = [0,0,0];
    delta_pos = [1e9; 1e9; 1e9; 0; 0];
    while ((norm(delta_pos(1:3,1))>1e-3))
        for idx_sv = 1 : length (sv_pos)
            pos1_pr0(idx_sv) = norm(sv_pos(idx_sv,2:4)-pos_ini);
            prn = sv_pos(idx_sv,1);
            if (prn > 32) && (prn <= 56)
                y(idx_sv,1) = pr(idx_sv) - pos1_pr0(idx_sv) - delta_pos(4,1);
                H(idx_sv,1:4) = [(pos_ini - sv_pos(idx_sv,2:4))./pos1_pr0(idx_sv),1];
            elseif prn > 86
                y(idx_sv,1) = pr(idx_sv) - pos1_pr0(idx_sv) - delta_pos(4,1) - delta_pos(5,1);
                H(idx_sv,1:5) = [(pos_ini - sv_pos(idx_sv,2:4))./pos1_pr0(idx_sv),1,1];
            end
        end
    delta_pos = inv(H'*H)*H'*y;
    pos_ini = pos_ini + delta_pos(1:3,1)';
    dt(1,2) = dt(1,2) + delta_pos(4,1);
    dt(1,3) = dt(1,3) + delta_pos(5,1);
    end
    pos_est = pos_ini;
    pr_resi = y - H*delta_pos;
    tag = 6;
end

function [pos_est,pr_resi] = LS_pos_NB(pr,sv_pos,pos_ini)

% Least square positioning (GPS& BDS)
% Input:    pr      =   pseudorange measurement
%           sv_pos  =   satellite position in ECEF
%           pos_ini =   initial guess of position in ECEF
%
% Output:   pos_est =   position estimation
%           sse     =   sum of square error

if size(pr,1)>=5
    dt = [0,0,0];
    delta_pos = [1e9; 1e9; 1e9; 0; 0];
    while ((norm(delta_pos(1:3,1))>1e-3))
        for idx_sv = 1 : length (sv_pos)
            pos1_pr0(idx_sv) = norm(sv_pos(idx_sv,2:4)-pos_ini);
            prn = sv_pos(idx_sv,1);
            if prn <= 32
                y(idx_sv,1) = pr(idx_sv) - pos1_pr0(idx_sv) - delta_pos(4,1);
                H(idx_sv,1:4) = [(pos_ini - sv_pos(idx_sv,2:4))./pos1_pr0(idx_sv),1];
            elseif prn > 86
                y(idx_sv,1) = pr(idx_sv) - pos1_pr0(idx_sv) - delta_pos(4,1) - delta_pos(5,1);
                H(idx_sv,1:5) = [(pos_ini - sv_pos(idx_sv,2:4))./pos1_pr0(idx_sv),1,1];
            end
        end
    delta_pos = inv(H'*H)*H'*y;
    pos_ini = pos_ini + delta_pos(1:3,1)';
    dt(1,1) = dt(1,1) + delta_pos(4,1);
    dt(1,3) = dt(1,3) + delta_pos(5,1);
    end
    pos_est = pos_ini;
    pr_resi = y - H*delta_pos;
    tag = 5;
end

function [pos_est,pr_resi] = LS_pos_NG(pr,sv_pos,pos_ini)

% Least square positioning (GPS & GLONASS only)
% Input:    pr      =   pseudorange measurement
%           sv_pos  =   satellite position in ECEF
%           pos_ini =   initial guess of position in ECEF
%
% Output:   pos_est =   position estimation
%           sse     =   sum of square error

if size(pr,1)>=5
    dt = [0,0,0];
    delta_pos = [1e9; 1e9; 1e9; 0; 0];
    while ((norm(delta_pos(1:3,1))>1e-3))
        for idx_sv = 1 : length (sv_pos)
            pos1_pr0(idx_sv) = norm(sv_pos(idx_sv,2:4)-pos_ini);
            prn = sv_pos(idx_sv,1);
            if prn <= 32
                y(idx_sv,1) = pr(idx_sv) - pos1_pr0(idx_sv) - delta_pos(4,1);
                H(idx_sv,1:4) = [(pos_ini - sv_pos(idx_sv,2:4))./pos1_pr0(idx_sv),1];
            elseif (prn > 32) && (prn <= 56)
                y(idx_sv,1) = pr(idx_sv) - pos1_pr0(idx_sv) - delta_pos(4,1) - delta_pos(5,1);
                H(idx_sv,1:5) = [(pos_ini - sv_pos(idx_sv,2:4))./pos1_pr0(idx_sv),1,1];
            end
        end
    delta_pos = inv(H'*H)*H'*y;
    pos_ini = pos_ini + delta_pos(1:3,1)';
    dt(1,1) = dt(1,1) + delta_pos(4,1);
    dt(1,2) = dt(1,2) + delta_pos(5,1);
    end
    pos_est = pos_ini;
    pr_resi = y - H*delta_pos;
    tag = 4;
end

function [pos_est,pr_resi] = LS_pos_NGB(pr,sv_pos,pos_ini)

% Least square positioning (GPS, GLONASS & BDS)
% Input:    pr      =   pseudorange measurement
%           sv_pos  =   satellite position in ECEF
%           pos_ini =   initial guess of position in ECEF
%
% Output:   pos_est =   position estimation
%           sse     =   sum of square error

if size(pr,1)>=6
    dt = [0,0,0];
    delta_pos = [1e9; 1e9; 1e9; 0; 0; 0];
    while ((norm(delta_pos(1:3,1))>1e-3))
        for idx_sv = 1 : length (sv_pos)
            pos1_pr0(idx_sv) = norm(sv_pos(idx_sv,2:4)-pos_ini);
            prn = sv_pos(idx_sv,1);
            if prn <= 32
                y(idx_sv,1) = pr(idx_sv) - pos1_pr0(idx_sv) - delta_pos(4,1);
                H(idx_sv,1:4) = [(pos_ini - sv_pos(idx_sv,2:4))./pos1_pr0(idx_sv),1];
            elseif (prn > 32) && (prn <= 56)
                y(idx_sv,1) = pr(idx_sv) - pos1_pr0(idx_sv) - delta_pos(4,1) - delta_pos(5,1);
                H(idx_sv,1:5) = [(pos_ini - sv_pos(idx_sv,2:4))./pos1_pr0(idx_sv),1,1];
            elseif (prn > 56) && (prn <= 123)
                y(idx_sv,1) = pr(idx_sv) - pos1_pr0(idx_sv) - delta_pos(4,1)  - delta_pos(6,1);
                H(idx_sv,1:6) = [(pos_ini - sv_pos(idx_sv,2:4))./pos1_pr0(idx_sv),1,0,1]; 
            end
        end
    delta_pos = inv(H'*H)*H'*y;
    pos_ini = pos_ini + delta_pos(1:3,1)';
    dt(1,1) = dt(1,1) + delta_pos(4,1);
    dt(1,2) = dt(1,2) + delta_pos(5,1);
    dt(1,3) = dt(1,3) + delta_pos(6,1);
    end
    pos_est = pos_ini;
    pr_resi = y - H*delta_pos;
    tag = 7;
end