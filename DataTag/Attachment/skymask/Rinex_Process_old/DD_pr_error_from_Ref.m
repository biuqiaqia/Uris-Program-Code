function [Pr_error] = DD_pr_error_from_Ref(filename_nav,ref_obs_path,ref_xyz,user_data,user_xyz)

% ---------------------------------------------------
% Function to calculate the pseudorange error from
% reference station .obs
% by GH.Zhang 2018/08/28
% ---------------------------------------------------

% clc;
% clear;
% close all;
% D2R = pi/180;
% R2D = 180/pi;
% %input
% filename_nav = 'data\eph\hkqt246j.18n';
% ref_obs_path = 'data\eph\hkqt246j.18o';
% ref_xyz = [-2421567.8916,5384910.5631,2404264.3943];
% user_llh = [22.314430,114.206308,5];
% user_xyz = llh2xyz([user_llh(1,1:2)*D2R,user_llh(1,3)]);
% 
% load('data\Airport_test\ubx_data.mat');
% user_data = data_all;
% clear data_all;
% user_xyz = ones(size(user_data,1),3).*user_xyz;
        
%% Processing Ref OBS file
D2R = pi/180;
R2D = 180/pi;
global weights 
weights = 0;
global FTABLE;
FTABLE=finv(0.9995,1,1:200)';

%enabled constellations
GPS_flag = 1;  GLO_flag = 1;  GAL_flag = 1;  BDS_flag = 1;  QZS_flag = 1; SBS_flag = 0;
GPS_col = 'b'; GLO_col = 'r'; GAL_col = 'g'; BDS_col = 'c'; QZS_col = 'm';

[constellations] = goGNSS.initConstellation(GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag);
nSatTot = constellations.nEnabledSat;

goGNSS;
                [Eph, iono] = load_RINEX_nav(filename_nav);            

                [pr1_R, ph1_R, pr2_R, ph2_R, dop1_R, dop2_R, snr1_R, snr2_R, ...
                 time_GPS, time_R, week_R, date_R, pos_R, interval, antoff_R, antmod_R] = ...
                 load_RINEX_obs(ref_obs_path, constellations);
             
fprintf('Ref RINEX files loading complete\n');

is_bias_tot=NaN(6,length(time_GPS));
[~, all_sow] = time2weektow(time_GPS);
    
% lon
idx_run_data = 0;
XR0 = [];
t=0;
id_ts = find(all_sow==user_data{1,1});
id_te = find(all_sow==user_data{end,1});


if isempty(id_te)
    id_te = find(all_sow==user_data{3600,1});
end



for t = id_ts:id_te
    XR=[];
    time_rx = time_GPS(t);       
    time_rx_ = mod(time_rx,604800);

    idx_run_data = idx_run_data + 1;
    pr1 = pr1_R(:,t);
    dop1 = dop1_R(:,t);  
    snr = snr1_R(:,t);
    nSatTot = size(pr1,1);
    Eph_t = rt_find_eph(Eph, time_GPS(t), nSatTot);      
    lambda = goGNSS.getGNSSWavelengths(Eph_t, nSatTot);
    eph_avail = Eph(30,:);
    dtR = 0;
    %% Rough Satellite Position Estimation
    [XS, dtS, XS_tx, VS_tx, time_tx, no_eph, sys] = satellite_positions(time_rx, 20200*1e3*ones(1,length(eph_avail)), eph_avail, Eph, [], [], zeros(1,length(eph_avail)), zeros(1,length(eph_avail)), dtR);

    sat = find(pr1 ~= 0);
    eph_avail = Eph(30,:);
    sat = sat(ismember(sat, eph_avail));
    
    rdata_all{t-id_ts+1,1} = time_rx_; % GPS TOW; 
    
        for jdx = 1 : length(XS)
                target_pos_xyz = [XS(jdx,1) XS(jdx,2) XS(jdx,3)];    
                target_pos_enu = xyz2enu(target_pos_xyz,ref_xyz);
                azimuth(jdx) = ((pi/2)-atan2(target_pos_enu(1),target_pos_enu(2)))*R2D;
                if (azimuth(jdx)<0)
                    azimuth(jdx) = azimuth(jdx) + 360 ;
                end
                elevation(jdx) = ( R2D*atan(target_pos_enu(3)/norm(target_pos_enu(1:2))));
        end
        rdata_all{t-id_ts+1,6} = eph_avail; % ALL PRN;
        rdata_all{t-id_ts+1,7} = elevation; % elevation angle for each grid location;
        rdata_all{t-id_ts+1,8} = azimuth; % azimuth angle for each grid location;                    
    
    %% WLS to estimate rough receiver position
    sat = find(pr1 ~= 0);
    eph_avail = Eph(30,:);
    sat = sat(ismember(sat, eph_avail));

    flag_gps_enough = 1;
    if (length(find(sat <= 32))<4)
        flag_gps_enough = 0;    
    end
    min_nsat = 4;
    
    XR_LLH(idx_run_data,:) = zeros(1,3);
    VR_LLH(idx_run_data,:) = zeros(1,3);
    Time_WLS(idx_run_data) = all_sow(t);
    
    if (size(sat,1) >= min_nsat && flag_gps_enough)
    
    lambda = lambda(sat);
    pseudorange = pr1(sat);
    doppler = dop1(sat);
    snr = snr(sat);
    
    [XS, dtS, XS_tx, VS_tx, time_tx, no_eph, sys] = satellite_positions(time_rx, pseudorange, sat, Eph, [], [], zeros(1,length(sat)), zeros(1,length(sat)), dtR);
       
    num_sys  = length(unique(sys(sys ~= 0)));
    min_nsat = 3 + num_sys;
    
    if (isempty(XR0))
        XR0 = zeros(3,1);
    end
    
    index = find(no_eph == 0);   
    nsat_avail = length(index);

    if (nsat_avail < min_nsat) %if available observations are not enough, return empty variables
        fprintf('not enough sv\n');
    else
        %iterative least-squares from XR0,i.e. given coordinates or the center of the Earth (i.e. [0; 0; 0])
        n_iter_max = 5;
        n_iter = 0;
        var_SPP(1) = Inf;
        SPP_threshold=4; %meters 
        while(var_SPP(1) > SPP_threshold^2 && n_iter < n_iter_max)
            [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num, bad_obs, bad_epoch, var_SPP] = LS_SA_code_Qmo(XR0, XS(index,:), pseudorange(index), snr(index), zeros(nsat_avail,1), zeros(nsat_avail,1), dtS(index), zeros(nsat_avail,1), zeros(nsat_avail,1), sys(index), SPP_threshold);       
            XR0 = XR;
            n_iter = n_iter + 1;
        end        
        %satellite topocentric coordinates (azimuth, elevation, distance)
        [az, el, dist_] = topocent(XR, XS(index,:));        
        %cartesian to geodetic conversion of ROVER coordinates
        [phiR, lamR, hR] = cart2geod(XR(1), XR(2), XR(3));
        %radians to degrees
        phiR = phiR * 180 / pi;
        lamR = lamR * 180 / pi;
        %computation of tropospheric errors
        err_tropo = tropo_error_correction(el, hR);
        %computation of ionospheric errors
        err_iono = iono_error_correction(phiR, lamR, az, el, time_rx, iono, []);  

        %accurate Clock estimation
        [dtR, var_dtR, bad_obs, bad_epoch, var_SPP, ~, is_bias] = LS_SA_code_clock_Qmo(pseudorange(index), snr(index), el, dist_, dtS(index), err_tropo, err_iono, sys(index), SPP_threshold);
        %accurate SV Pos
        [XS, dtS, XS_tx, VS_tx, time_tx, no_eph, sys] = satellite_positions(time_rx, pseudorange, sat, Eph, [], [], zeros(1,length(sat)), zeros(1,length(sat)), dtR(1));
        
        %SV velocity estimation (learned from RTKLIB)
        tt = 1e-3;
        [XS_tt, dtS_tt, XS_tx_tt, VS_tx_tt, time_tx_tt, no_eph, sys] = satellite_positions(time_rx+tt, pseudorange, sat, Eph, [], [], zeros(1,length(sat)), zeros(1,length(sat)), dtR(1));        
        VS = (XS_tt - XS) / tt;
        dtSV = (dtS_tt - dtS) / tt;
        
        %aacurate rover pos estimation
        n_iter = 0;
        outlier_prn{t,1} = [];
        while(var_SPP(1) > SPP_threshold^2 && n_iter < n_iter_max)
            [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num, bad_obs, bad_epoch, var_SPP] = LS_SA_code_Qmo(XR0, XS(index,:), pseudorange(index), snr(index), el, zeros(nsat_avail,1), dtS(index), err_tropo, err_iono, sys(index), SPP_threshold);       
            outlier_prn{t,1} = bad_obs;
            XR0 = XR;
            n_iter = n_iter + 1;
            outlier_id{t,n_iter} = bad_epoch;
        end
        %satellite topocentric coordinates (azimuth, elevation, distance)
        [az, el, dist_] = topocent(XR, XS(index,:)); 
        az_enu = az.*(-1)+90;
        az_enu(az_enu<0,1) = az_enu(az_enu<0,1)+360;
        
        %aacurate rover vel estimation
%         [VR, dtRV, pr_rate_lt,A,LOS_vector_lt] = LS_SA_code_Vel_Qmo_debug(XR, XS(index,:), VS(index,:), doppler(index), snr(index), el, lambda(index), sys(index));       

    end    
    end
    
    if ~isempty(XR)
    rcr_loc_temp = xyz2llh(XR);
    rdata_all{t-id_ts+1,2} = [rcr_loc_temp(2)*R2D rcr_loc_temp(1)*R2D rcr_loc_temp(3)]; % location that receiver estimated (LLA)
    rdata_all{t-id_ts+1,3} = sat(index,1); % PRN list that receiver received 
    rdata_all{t-id_ts+1,4} = el; % elevation angle that receiver received
    rdata_all{t-id_ts+1,5} = [az_enu,az]; % azimuth angle that receiver received [ENU,NED]
    rdata_all{t-id_ts+1,9} = snr(index); % Carrier to noise ratio that receiver received 
    rdata_all{t-id_ts+1,10} = [sat(index,1),pseudorange(index,1)+dtS(index,1).*goGNSS.V_LIGHT-err_tropo-err_iono,XS(index,:)];% Measurement data
    rdata_all{t-id_ts+1,11} = [doppler(index,1),VS(index,:),lambda(index,1)];% doppler shift
    end
    % The index here excludes some satellite without ephemeris
end
fprintf('Ref RINEX data process complete\n');

%% DD Pseudorange Error
for idu = 1:1:min(size(user_data,1),3600)
    gps_s = user_data{idu,1};
    idr = find([rdata_all{:,1}]==gps_s);
    if ~isempty(rdata_all{idr,10})
    ref_sdata{idu,1} = [];
    user_sdata{idu,1} = [];
    for upr_id = 1:1:size(user_data{idu,10},1)
        rpr_id = find(rdata_all{idr,10}(:,1)==user_data{idu,10}(upr_id,1));
        if ~isempty(rpr_id)
            user_sdata{idu,1} = [user_sdata{idu,1};user_data{idu,10}(upr_id,:),user_data{idu,4}(upr_id,1)];
            ref_sdata{idu,1} = [ref_sdata{idu,1};rdata_all{idr,10}(rpr_id,:),rdata_all{idr,4}(rpr_id,1)];
        end
    end
    [~,y_dd,H_dd,~,masterSV] = DD_pr_error_cal(ref_sdata{idu,1},user_sdata{idu,1},ref_xyz);
    delta_pos = ref_xyz-user_xyz(idu,1:3);
    if ~isempty(H_dd)
    D = H_dd*delta_pos';
    delta_Y = D-y_dd;
    Pr_error{idu,1} = [];
    idy = 1;
    for idm = 1:1:size(ref_sdata{idu,1},1)
        if ref_sdata{idu,1}(idm,1)>=56 && ref_sdata{idu,1}(idm,1)<=86 % not support galileo
            Pr_error{idu,1}(idm,1) = nan;
        else
            if ~isempty(masterSV{1,1})&&idm == masterSV{1,1}
                Pr_error{idu,1}(idm,1) = 0;
            elseif ~isempty(masterSV{1,2})&&idm == masterSV{1,2}
                Pr_error{idu,1}(idm,1) = 0;
            elseif ~isempty(masterSV{1,3})&&idm == masterSV{1,3}
                Pr_error{idu,1}(idm,1) = 0;
            else
                Pr_error{idu,1}(idm,1) = delta_Y(idy);
                idy = idy + 1;
            end
        end
    end
    Pr_error{idu,1}=[ref_sdata{idu,1}(:,1),Pr_error{idu,1},ref_sdata{idu,1}(:,6)];
    Pr_error{idu,2}=masterSV;
    end
    else
        Pr_error{idu,1} = [user_data{idu,10}(:,1),nan(size(user_data{idu,10},1),2)];
    end
end





















































