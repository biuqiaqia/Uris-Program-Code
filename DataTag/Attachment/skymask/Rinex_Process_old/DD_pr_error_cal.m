function [delta_pos_dd,y_dd,H_dd,DD_DOP,masterSV] = DD_pr_error_cal(ref_sdata,oth_sdata,ref_pos)

% Least square positioning (conventional) with elevation selection
% Input:    ref_sdata   =   ego vehicle satellite data
%           oth_sdata   =   another vehicle satellite data
%
% Output:   delta_pos_dd    =   distance between ego and another vehicle
%
% copy from LT

if size(ref_sdata,1) > 2
    idn = [];
    idg = [];
    idb = [];
    for ids = 1:1:size(ref_sdata,1)
        if ref_sdata(ids,1) <= 32
            idn = [idn;ids];
        elseif ref_sdata(ids,1) > 32 && ref_sdata(ids,1) <= 56
            idg = [idg;ids];
        elseif ref_sdata(ids,1) > 86
            idb = [idb;ids];
        end
    end
    
idx_sv_dd = 0;

refv_spos = ref_sdata(:,3:5);% ego vehicle satelite positions
pos1_pr = ref_sdata(:,2);% ego vehicle pseudoranges
pos2_pr = oth_sdata(:,2);% other vehicle pseudoranges
refPos_DD = ref_pos;% estimated location of reference
sum_ele = ref_sdata(:,6)+oth_sdata(:,6);

if ~isempty(idn)% GPS DD
    idn_masterSV = find(sum_ele(:,1)==max(sum_ele(idn,1)));
    pos1_prm = norm(ref_sdata(idn_masterSV,3:5)-refPos_DD);
    for idx_sv = idn'
        if idx_sv ~= idn_masterSV
            idx_sv_dd = idx_sv_dd + 1;
            pos1_pr0 = norm(ref_sdata(idx_sv,3:5)-refPos_DD);
            y_dd(idx_sv_dd,1) = (pos1_pr(idx_sv) - pos2_pr(idx_sv)) - (pos1_pr(idn_masterSV) - pos2_pr(idn_masterSV));
            H_dd(idx_sv_dd,:) = (refPos_DD - refv_spos(idx_sv,:))./pos1_pr0 - (refPos_DD - refv_spos(idn_masterSV,:))./pos1_prm;
            ref_sv_enu = xyz2enu(refv_spos(idn_masterSV,:),refPos_DD);
            oth_sv_enu = xyz2enu(refv_spos(idx_sv,:),refPos_DD);
            A(idx_sv_dd,1:3) = [oth_sv_enu'./pos1_pr0-ref_sv_enu'./pos1_prm];
        end
    end
else
    idn_masterSV = [];
end


if ~isempty(idg)% GLONASS DD
    idg_masterSV = find(sum_ele(:,1)==max(sum_ele(idg,1)));
    pos1_prm = norm(ref_sdata(idg_masterSV,3:5)-refPos_DD);
    for idx_sv = idg'
        if idx_sv ~= idg_masterSV
            idx_sv_dd = idx_sv_dd + 1;
            pos1_pr0 = norm(ref_sdata(idx_sv,3:5)-refPos_DD);
            y_dd(idx_sv_dd,1) = (pos1_pr(idx_sv) - pos2_pr(idx_sv)) - (pos1_pr(idg_masterSV) - pos2_pr(idg_masterSV));
            H_dd(idx_sv_dd,:) = (refPos_DD - refv_spos(idx_sv,:))./pos1_pr0 - (refPos_DD - refv_spos(idg_masterSV,:))./pos1_prm;
            ref_sv_enu = xyz2enu(refv_spos(idg_masterSV,:),refPos_DD);
            oth_sv_enu = xyz2enu(refv_spos(idx_sv,:),refPos_DD);
            A(idx_sv_dd,1:3) = [oth_sv_enu'./pos1_pr0-ref_sv_enu'./pos1_prm];
        end
    end
else
    idg_masterSV = [];
end

if ~isempty(idb)% BDS DD
    idb_masterSV = find(sum_ele(:,1)==max(sum_ele(idb,1)));
    pos1_prm = norm(ref_sdata(idb_masterSV,3:5)-refPos_DD);
    for idx_sv = idb'
        if idx_sv ~= idb_masterSV
            idx_sv_dd = idx_sv_dd + 1;
            pos1_pr0 = norm(ref_sdata(idx_sv,3:5)-refPos_DD);
            y_dd(idx_sv_dd,1) = (pos1_pr(idx_sv) - pos2_pr(idx_sv)) - (pos1_pr(idb_masterSV) - pos2_pr(idb_masterSV));
            H_dd(idx_sv_dd,:) = (refPos_DD - refv_spos(idx_sv,:))./pos1_pr0 - (refPos_DD - refv_spos(idb_masterSV,:))./pos1_prm;
            ref_sv_enu = xyz2enu(refv_spos(idb_masterSV,:),refPos_DD);
            oth_sv_enu = xyz2enu(refv_spos(idx_sv,:),refPos_DD);
            A(idx_sv_dd,1:3) = [oth_sv_enu'./pos1_pr0-ref_sv_enu'./pos1_prm];
        end
    end
else
    idb_masterSV = [];
end



% for idx_sv = 1 : size (refv_spos,1)
%     if idx_sv ~= idx_masterSV
%         idx_sv_dd = idx_sv_dd + 1;
% 
%         y_dd(idx_sv_dd,1) = (pos1_pr(idx_sv) - pos2_pr(idx_sv)) - (pos1_pr(idx_masterSV) - pos2_pr(idx_masterSV));
%         H_dd(idx_sv_dd,:) = (refPos_DD - refv_spos(idx_sv,:))./pos1_pr(idx_sv) - (refPos_DD - refv_spos(idx_masterSV,:))./pos1_pr(idx_masterSV);
%     end
% end
    if idx_sv_dd > 2
        masterSV{1} = idn_masterSV;
        masterSV{2} = idg_masterSV;
        masterSV{3} = idb_masterSV;        
        
        delta_pos_dd = H_dd\y_dd;
        delta_pos_dd = delta_pos_dd';
        Q = (A'*A)^-1;
        DD_DOP(1,1) = sqrt(Q(1,1)+Q(2,2));
        DD_DOP(1,2) = sqrt(Q(3,3));
    else
        masterSV = [];
        delta_pos_dd = [];
        DD_DOP = [0,0];
    end
else
    masterSV = [];
    y_dd=[];
    H_dd=[];
    delta_pos_dd = [];
    DD_DOP = [0,0];
end