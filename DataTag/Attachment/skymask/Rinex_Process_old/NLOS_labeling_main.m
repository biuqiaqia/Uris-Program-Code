% ---------------------------------------------------
% Main code for NLOS labeling
% based on ground truth skymask
% features are included
% Updated with DD Pr_Error
% Corrected with NED skymask and ephemeris azimuth
% by GH.Zhang 2018/09/14
% ---------------------------------------------------

clc;
clear all;
close all;

D2R = pi/180;
R2D = 180/pi;
m2lat = 1/110734;
m2lon = 1/103043;

%% File Input
% OPPO data
filename_obs = 'COM6___9600_230709_084159';% Ublox
% load(['data\',filename_obs,'_nmea.mat']);
% gt_waypoint = [rcvPosResult(1,2:4),1;...
%                rcvPosResult(end,2:4),0];%LLH and epoch
gt_waypoint = [22.2977444,114.179022,3,1;...
               22.2977444001,114.179022001,3.001,5000];
filename_nav = 'data\hksc190i.23n';           
ref_obs_path = 'data\hksc190i.23o';



% Processing Setup        
% file_path = 'data\building_model\demo_HK_TST_fix.kml';
% skymask_path = 'data\building_model\skymask_TST_ned';
ref_xyz = [-2414266.9197,5386768.9868,2407460.0314];%hksc
existing_skymask_tag = 1; % 0-no/ 1-exist/ 2-all open
skyplot_output = 0;
DD_pr_error_tag = 2;
% GNSS_type == 1;
NLOS_labeling_tag = 2;
        
%% Read Eph and Obs data
obs_path = ['data\',filename_obs,'.obs'];
global weights 
weights = 0;
global FTABLE;
FTABLE=finv(0.9995,1,1:200)';

%enabled constellations
GPS_flag = 1;  GLO_flag = 1;  GAL_flag = 1;  BDS_flag = 1;  QZS_flag = 1; SBS_flag = 0;
GPS_col = 'b'; GLO_col = 'r'; GAL_col = 'g'; BDS_col = 'c'; QZS_col = 'm';


[constellations] = goGNSS.initConstellation(GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag);
% constellations.nEnabledSat = constellations.nEnabledSat + 40;
nSatTot = constellations.nEnabledSat;

goGNSS;
                [Eph, iono] = load_RINEX_nav(filename_nav);            

                [pr1_R, ph1_R, pr2_R, ph2_R, dop1_R, dop2_R, snr1_R, snr2_R, ...
                 time_GPS, time_R, week_R, date_R, pos_R, interval, antoff_R, antmod_R] = ...
                 load_RINEX_obs(obs_path, constellations);
 

pr1_R(128:end,:)=[];
pr1_R([90,92],:) = 0;

fprintf('RINEX files loading complete\n');

%% Main start
user_pos_xyz0 = NED_to_ECEF_pos(gt_waypoint(1,1),gt_waypoint(1,2),gt_waypoint(1,3));
is_bias_tot=NaN(6,length(time_GPS));
[~, all_sow] = time2weektow(time_GPS);
disp(['---> Epoch Amount: ',num2str(size(all_sow,1))]);
    
% lon
idx_run_data = 0;
XR0 = [];
t=0;
id_data_1=1;

for t = 1:length(time_GPS)
    XR=[];
%     t = t + 1;
    time_rx = time_GPS(t);       
    time_rx_ = mod(time_rx,604800);

    idx_run_data = idx_run_data + 1;    
    if (idx_run_data == 1)
%         fprintf('Processing...\n');
    end      
    
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
    
    %---massive pr exclusion (HUAWEI)---
%     id_mpr = [];
%     for idpr = 1:1:size(pr1,1)
%         if (pr1(idpr,1) > 4e7 || pr1(idpr,1) < 1.5e7) && pr1(idpr,1) ~= 0 
%             id_mpr = [id_mpr;idpr];
%         else
%             
%         end
%     end
%     prn_bug{t,1} = id_mpr;
%     pr1(id_mpr,1)=0;
    %---end exclusion---

    
    sat = find(pr1 ~= 0);
    eph_avail = Eph(30,:);
    sat = sat(ismember(sat, eph_avail));
    
    data_all{t,1} = time_rx_; % GPS TOW; 
    
        for jdx = 1 : length(XS)
                target_pos_xyz = [XS(jdx,1) XS(jdx,2) XS(jdx,3)];    
                target_pos_enu = xyz2enu(target_pos_xyz,user_pos_xyz0);
                azimuth(jdx) = ((pi/2)-atan2(target_pos_enu(1),target_pos_enu(2)))*R2D;
                if (azimuth(jdx)<0)
                    azimuth(jdx) = azimuth(jdx) + 360 ;
                end
                elevation(jdx) = ( R2D*atan(target_pos_enu(3)/norm(target_pos_enu(1:2))));
        end
        eph_az_ned = azimuth.*(-1)+90;
        eph_az_ned(1,eph_az_ned<0) = eph_az_ned(1,eph_az_ned<0)+360;
        
        
        data_all{t,6} = eph_avail; % ALL PRN;
        data_all{t,7} = elevation; % elevation angle for each grid location;
        data_all{t,8} = eph_az_ned; % azimuth angle for each grid location;                    
    
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
        [VR, dtRV, pr_rate_lt,A,LOS_vector_lt] = LS_SA_code_Vel_Qmo_debug(XR, XS(index,:), VS(index,:), doppler(index), snr(index), el, lambda(index), sys(index));       

    end    
    end
    
    if ~isempty(XR)
    rcr_loc_temp = xyz2llh(XR);
    data_all{t,2} = [rcr_loc_temp(2)*R2D rcr_loc_temp(1)*R2D rcr_loc_temp(3)]; % location that receiver estimated (LLA)
    data_all{t,3} = sat(index,1); % PRN list that receiver received 
    data_all{t,4} = el; % elevation angle that receiver received
    data_all{t,5} = [az_enu,az]; % azimuth angle that receiver received [ENU,NED]
    data_all{t,9} = snr(index); % Carrier to noise ratio that receiver received 
    data_all{t,10} = [sat(index,1),pseudorange(index,1)+dtS(index,1).*goGNSS.V_LIGHT-err_tropo-err_iono,XS(index,:)];% Measurement data
    data_all{t,11} = [doppler(index,1),VS(index,:),lambda(index,1)];% doppler shift
    data_all{t,12} = VR';
    end
    % The index here excludes some satellite without ephemeris
    % save only 1Hz data
    if rem(all_sow(t),1)==0
        for idc = 1:1:11
        data_all_1{id_data_1,idc}=data_all{t,idc};
        end
        id_data_1=id_data_1+1;
    end
    % Available SV & LOS Ratio
    [~,ia,~] = unique(data_all{t,6}); % negelct the duplicated ephs
    idx_eph_visible_sat = find((data_all{t,6}(ia)>=87 | (data_all{t,6}(ia)<=56)) & data_all{t,7}(ia)>0);
    data_all{t,6} = data_all{t,6}(ia(idx_eph_visible_sat));
    data_all{t,7} = data_all{t,7}(ia(idx_eph_visible_sat));
    data_all{t,8} = data_all{t,8}(ia(idx_eph_visible_sat));
end

% Delete Novatel 5Hz extra data
% outlier_prn(rem(all_sow,1)~=0,:)=[];
% all_sow(rem(all_sow,1)~=0,:)=[];
% data_all = data_all_1;
% disp(['---> Epoch Amount Refine to: ',num2str(size(all_sow,1))]);
% t = size(all_sow,1);


% clearvars -except all_sow data_all D2R R2D m2lat m2lon gt_waypoint ...
%                   ref_obs_path ref_xyz file_path skymask_path t ...
%                   existing_skymask_tag skyplot_output DD_pr_error_tag ...
%                   outlier_prn user_pos_xyz0 filename_nav filename_obs;


%% Ground truth
[num_data_epoch, ~] = size(data_all);
% gt_waypoint
gt_lon_epoch = [];
gt_lat_epoch = [];
gt_waypoint(end,4)=t;
for idwp = 1:1:size(gt_waypoint,1)-1
    gt_lon_epoch = [gt_lon_epoch,gt_waypoint(idwp,2):(gt_waypoint(idwp+1,2)-gt_waypoint(idwp,2))/(gt_waypoint(idwp+1,4)-gt_waypoint(idwp,4)):gt_waypoint(idwp+1,2)];
    gt_lat_epoch = [gt_lat_epoch,gt_waypoint(idwp,1):(gt_waypoint(idwp+1,1)-gt_waypoint(idwp,1))/(gt_waypoint(idwp+1,4)-gt_waypoint(idwp,4)):gt_waypoint(idwp+1,1)];
    gt_lon_epoch(end) = [];
    gt_lat_epoch(end) = [];
end
gt_lon_epoch = [gt_lon_epoch,gt_waypoint(end,2)];
gt_lat_epoch = [gt_lat_epoch,gt_waypoint(end,1)];

for id_gt = 1:1:num_data_epoch
    gt_xyz(id_gt,:) = NED_to_ECEF_pos(gt_lat_epoch(id_gt),gt_lon_epoch(id_gt),gt_waypoint(1,3));
end

%% Positioning Error from GT
for idt = 1:1:size(data_all,1)
    if ~isempty(data_all{idt,2})
    wls_llh(idt,1) = data_all{idt,2}(2);
    wls_llh(idt,2) = data_all{idt,2}(1);
    wls_llh(idt,3) = data_all{idt,2}(3);
    wls_xyz(idt,:) = NED_to_ECEF_pos(wls_llh(idt,1),wls_llh(idt,2),gt_waypoint(1,3));
    Error.WLS_2D(idt,1) = norm(gt_xyz(idt,:)-wls_xyz(idt,:));
    wls_xyz(idt,:) = NED_to_ECEF_pos(wls_llh(idt,1),wls_llh(idt,2),wls_llh(1,3));
    Error.WLS_3D(idt,1) = norm(gt_xyz(idt,:)-wls_xyz(idt,:));
    % Pseudorange_residual
    if ~isempty(outlier_prn{idt,1})&&1
        data_all{idt,3}(outlier_prn{idt,1},:)=[];
        data_all{idt,4}(outlier_prn{idt,1},:)=[];
        data_all{idt,5}(outlier_prn{idt,1},:)=[];
        data_all{idt,9}(outlier_prn{idt,1},:)=[];
        data_all{idt,10}(outlier_prn{idt,1},:)=[];
        data_all{idt,11}(outlier_prn{idt,1},:)=[];
    end
    [ls_xyz(idt,:),pr_residual{idt,1},tag] = LS_mpos_pr_resi(data_all{idt,10}(:,2),data_all{idt,10}(:,[1,3:5]),user_pos_xyz0);
    ls_llh(idt,:) = xyz2llh(ls_xyz(idt,:));
    ls_llh(idt,1:2) = ls_llh(idt,1:2).*R2D;
    pr_residual{idt,2} = pr_residual{idt,1}./mean(abs(pr_residual{idt,1}));
    % LOS_ratio
    [unique_prn,ia,~] = unique(data_all{idt,6}); % negelct the duplicated ephs        
    idx_eph_visible_sat = find((data_all{idt,6}(ia)>=87 | (data_all{idt,6}(ia)<=56)) & data_all{idt,7}(ia)>0);
    temp_prn_eph = data_all{idt,6}(ia(idx_eph_visible_sat));
    LOS_ratio(idt,1) = size(data_all{idt,4},1)/size(temp_prn_eph,2);
    end
end
disp(['2D Error -> ',num2str(mean(Error.WLS_2D)),'; 3D Error -> ',num2str(mean(Error.WLS_3D))]);

%% Single difference labelling
% for idt = 1:size(data_all,1)
%     n_id = data_all{idt,3}<=32;
%     [~,n_m_id] = max(data_all{idt,4}(n_id,:));
%     g_id = data_all{idt,3}>32 & data_all{idt,3}<=56;
%     [~,n_m_id] = max(data_all{idt,4}(n_id,:));
% end
%% NLOS labeling with skymask
% if existing_skymask_tag == 1
%     load(skymask_path);
%     for idt = 1:1:size(data_all,1)
%         disp(['---Processing Epoch ',num2str(idt),'/',num2str(size(all_sow,1)),'---']);
%         if ~isempty(data_all{idt,3})
%         [min_dlat idr_lat]= min(abs(lat_grid(1,:)-gt_lat_epoch(idt)));
%         [min_dlon idr_lon]= min(abs(lon_grid(:,1)-gt_lon_epoch(idt)));
%         gt_skymask{idt,1} = skymask{idr_lon,idr_lat}';
%         SV_elaz{idt,1} = [data_all{idt,3},data_all{idt,4},data_all{idt,5}(:,2)];
%         Label{idt,1} = ones(size(SV_elaz{idt,1},1),1);
%         for idm = 1:1:size(SV_elaz{idt,1},1)
%     %         SV_elaz(idm,4)=skymask(round(SV_elaz(idm,3))+1,2);
%             if gt_skymask{idt,1}(1,round(SV_elaz{idt,1}(idm,3))+1)>SV_elaz{idt,1}(idm,2)
%                 Label{idt,1}(idm,1)=0;
%             end
%         end
%         SV_elaz{idt,1}(:,4)=Label{idt,1};
%         end
%     end
% elseif existing_skymask_tag == -7
%     for idt = 1:1:size(data_all,1)
%         disp(['---Processing Epoch ',num2str(idt),'/',num2str(size(all_sow,1)),'---']);
%         temp_skymask = funtion3_map360(gt_waypoint(1,2),gt_waypoint(1,1),gt_waypoint(1,3),file_path,0);
%         gt_skymask{idt,1} = temp_skymask(:,2)';
%         SV_elaz{idt,1} = [data_all{idt,3},data_all{idt,4},data_all{idt,5}(:,2)];
%         Label{idt,1} = ones(size(SV_elaz{idt,1},1),1);
%         for idm = 1:1:size(SV_elaz{idt,1},1)
%     %         SV_elaz(idm,4)=skymask(round(SV_elaz(idm,3))+1,2);
%             if gt_skymask{idt,1}(1,round(SV_elaz{idt,1}(idm,3))+1)>SV_elaz{idt,1}(idm,2)
%                 Label{idt,1}(idm,1)=0;
%             end
%         end
%         SV_elaz{idt,1}(:,4)=Label{idt,1};
%     end
% elseif existing_skymask_tag == 2
%     for idt = 1:1:size(data_all,1)
%         disp(['---Processing Epoch ',num2str(idt),'/',num2str(size(all_sow,1)),'---']);
% %         temp_skymask = funtion3_map360(gt_waypoint(1,2),gt_waypoint(1,1),gt_waypoint(1,3),file_path,0);
%         gt_skymask{idt,1} = zeros(361,1)';
%         SV_elaz{idt,1} = [data_all{idt,3},data_all{idt,4},data_all{idt,5}(:,2)];
%         Label{idt,1} = ones(size(SV_elaz{idt,1},1),1);
%         for idm = 1:1:size(SV_elaz{idt,1},1)
%     %         SV_elaz(idm,4)=skymask(round(SV_elaz(idm,3))+1,2);
%             if gt_skymask{idt,1}(1,round(SV_elaz{idt,1}(idm,3))+1)>SV_elaz{idt,1}(idm,2)
%                 Label{idt,1}(idm,1)=0;
%             end
%         end
%         SV_elaz{idt,1}(:,4)=Label{idt,1};
%     end
% end

%% Pseudorange Rate Consistency
% From Doppler
% Pr_rate{1,1} = [ones(size(data_all{1,11},1),1).*9999,zeros(size(data_all{1,11},1),1)];
% for idt = 2:1:size(data_all,1)
%     for idm = 1:1:size(data_all{idt,11},1)
%         if ~isempty(data_all{idt,3})&&~isempty(data_all{idt-1,3})&&any(data_all{idt,3}(idm,1)==data_all{idt-1,3}(:,1))
%             temp_R = norm(data_all{idt,10}(idm,3:5)-wls_xyz(idt,:));
%             temp_LOS_Vector = (wls_xyz(idt,:) - data_all{idt,10}(idm,3:5))./temp_R;
%             Pr_rate{idt,1}(idm,1) = data_all{idt,11}(idm,1)*data_all{idt,11}(idm,5);            
%             id_pr_last = find(data_all{idt,10}(idm,1)==data_all{idt-1,10}(:,1));
%             Pr_rate{idt,1}(idm,2) = -(data_all{idt,10}(idm,2) - data_all{idt-1,10}(id_pr_last,2));
%         else
%             Pr_rate{idt,1}(idm,1) = 9999;
%             Pr_rate{idt,1}(idm,2) = 0;
%         end
%     end
%     Pr_rate{idt,1}(:,3) = Pr_rate{idt,1}(:,1)-Pr_rate{idt,1}(:,2);
% end

%% Pseudorange Error from Double Difference
if DD_pr_error_tag == 1
    disp('-----calculating Pseudorange Error from DD-----');
    [Pr_error] = DD_pr_error_from_Ref(filename_nav,ref_obs_path,ref_xyz,data_all,gt_xyz);
end

%% Labeling CSV Output
% disp('-----Writing CSV-----');
% out_dir = 'data\labels\';
% fid_out = fopen([out_dir,filename_obs,'_label.csv'],'w+');
% fprintf(fid_out,'GPS_Time(s),PRN,nSV,LOS_ratio,C/N0,Elevation,Azimuth,Pseudorange_residual,Normalized_Pseudorange_residual,Pr_rate_consitency,LOS/NLOS_label\n');
% for idt = 1:1:size(data_all,1)
%     for idm = 1:1:size(data_all{idt,3},1)
%         fprintf(fid_out,'%6d,%3d,%2d,%5.3f,%3.1f,%5.2f,%6.2f,%9.3f,%9.3f,%12.4f,%1d',...
%                     data_all{idt,1},...%GPS_Time
%                     data_all{idt,3}(idm),...%PRN
%                     size(data_all{idt,3},1),...%nSV
%                     LOS_ratio(idt,1),...%LOS_ratio
%                     data_all{idt,9}(idm),...%C/N0
%                     data_all{idt,4}(idm),...%Elevation
%                     data_all{idt,5}(idm,2),...%Azimuth
%                     pr_residual{idt,1}(idm),...%Pseudorange_residual
%                     pr_residual{idt,2}(idm),...%Normalized Pseudorange_residual
%                     Pr_rate{idt,1}(idm,1)-Pr_rate{idt,1}(idm,2),...%Delta_Pseudorange_Rate
%                     Label{idt,1}(idm)...%LOS/NLOS label
%                     );
%         if DD_pr_error_tag == 1
%             ide = find(Pr_error{idt,1}(:,1)==data_all{idt,3}(idm));
%             if ~isempty(ide)
%                 fprintf(fid_out,',%9.3f\n',Pr_error{idt,1}(ide,2));
%             else
%                 fprintf(fid_out,',%9.3f\n',999);
%             end
%         else
%             fprintf(fid_out,'\n');
%         end
%     end
% end
% fclose(fid_out);



%% NLOS Error Labeling (currently static)

for idt = 1:1:size(data_all,1)
    %         disp(['Skymasking...',num2str(idt),'/',num2str(size(data_all,1))]);
    %         [skymask{idt,1},~,~] = generate_skymask(gt_llh(idt,:), building_kml, [], 0.01);
    if ~isempty(data_all{idt,3}) && isequal(data_all{idt,3},Pr_error{idt,1}(:,1))
        %PrError
        data_all{idt,9}(:,5) = Pr_error{idt,1}(:,2);
        %skymask_vis
        data_all{idt,9}(:,6) = data_all{idt,4}>gt_skymask{idt,1}(ceil(data_all{idt,5}(:,2)))';
    else
        if isempty(data_all{idt,3})
            disp(['idt = ',num2str(idt),' empty!']);
            data_all{idt,9}(:,5:6) = nan(size(data_all{idt,9},1),2);
        else
            disp(['idt = ',num2str(idt),' not match!!!']);
            %PrError
            data_all{idt,9}(~ismember(data_all{idt,3},Pr_error{idt,1}(:,1)),5) = nan;
            data_all{idt,9}(ismember(data_all{idt,3},Pr_error{idt,1}(:,1)),5) = Pr_error{idt,1}(:,2);
            %skymask_vis
            data_all{idt,9}(:,6) = data_all{idt,4}>gt_skymask{idt,1}(ceil(data_all{idt,5}(:,2)))';
        end
        
    end
end

% save labeled measurements
cdata = data_all(:,1:8);
for idt = 1:1:size(data_all,1)
    for idm = 1:1:size(data_all{idt,3},1)
        cdata{idt,9}(idm,1) = data_all{idt,9}(idm,1);%C/N0
        cdata{idt,9}(idm,2) = pr_residual{idt,1}(idm);%Pseudorange_residual
        cdata{idt,9}(idm,3) = pr_residual{idt,2}(idm);%Normalized Pseudorange_residual
%         cdata{idt,9}(idm,4) = Pr_rate{idt,1}(idm,1)-Pr_rate{idt,1}(idm,2);%Delta_Pseudorange_Rate
        cdata{idt,9}(idm,5:6) = data_all{idt,9}(idm,5:6);%C/N0
%         if DD_pr_error_tag == 1
%             ide = find(Pr_error{idt,1}(:,1)==data_all{idt,3}(idm));
%             if ~isempty(ide)
%                 cdata{idt,9}(idm,5) = Pr_error{idt,1}(ide,2);% Pseudorange Error
%             else
%                 cdata{idt,9}(idm,5) = nan;
%             end
%         end
%         cdata{idt,9}(idm,6) = Label{idt,1}(idm);% Label
    end
end
cdata(:,10) = data_all(:,10);
cdata(:,11) = data_all(:,12);
gt_epoch = [gt_lat_epoch',gt_lon_epoch'];

%% Epoch selection
id_st = find([cdata{:,1}]==epoch_st);
id_ed = find([cdata{:,1}]==epoch_end);
cdata = cdata(id_st:id_ed,:);
gt_epoch = gt_epoch(id_st:id_ed,:);
ls_llh = ls_llh(id_st:id_ed,:);
Error.WLS_2D = Error.WLS_2D(id_st:id_ed,:);  
Error.WLS_3D = Error.WLS_3D(id_st:id_ed,:);  

save(['data\mat\',filename_obs,'.mat'],'cdata','gt_epoch');

% for idt = 1:1:size(data_all,1)
%     labeled_mea{idt,1} = data_all{idt,1};
%     labeled_mea{idt,2} = [data_all{idt,10},Label{idt,1}];
% %     labeled_mea{idt,3} = Label{idt,1};
% end
% save(['result\labeled_measurements\',filename_obs,'_mea.mat'],'labeled_mea','gt_xyz');





%% Plot
% WLS
ls_llh(ls_llh(:,1)==0,:)=nan;
% wls_llh(wls_llh(:,1)==0,:)=nan;

figure(1);
hold on;
plot(ls_llh(:,2),ls_llh(:,1),'b.','MarkerSize',25);
% plot(wls_llh(:,2),wls_llh(:,1),'r.','MarkerSize',25);
% plot(gt_epoch(:,2),gt_epoch(:,1),'g.-','MarkerSize',25);
hBase = plot_openstreetmap('Alpha',1,'Scale',50000,'BaseUrl', "http://a.tile.openstreetmap.org");
plot(ls_llh(:,2),ls_llh(:,1),'b.','MarkerSize',25);
% plot(wls_llh(:,2),wls_llh(:,1),'r.','MarkerSize',25);
% plot(gt_epoch(:,2),gt_epoch(:,1),'g.-','MarkerSize',25);

%Skyplot
if skyplot_output == 1
for idt = 1:1:size(data_all,1)
    skymask_ref = gt_skymask{idt,1};

    ref_data = [];
    for az = 1:1:360
        mask_el = skymask_ref(1,az);
        for el = 0:1:90
            if el <= mask_el
                tag = 2;
            else
                tag = 1;
            end
            ref_data = [ref_data;el,az,tag];
        end
    end

    elevation_ref = ref_data(:,1);
    azimuth_ref = ref_data(:,2);
    signaltype_ref = ref_data(:,3);
    
    h_gps = figure(idt+100);
    a = 1;
    patch([-a -a a  a],...
          [-a  a a -a],[239/255 239/255 239/255],'EdgeColor','none');

    radius = 0.9;
    LOSidx = find(signaltype_ref==1);
    NLOSidx = find(signaltype_ref==2 | signaltype_ref==3);
    hold on;
    h=plot(radius*(ones(size(LOSidx,1),1)-elevation_ref(LOSidx)./90).*sind(azimuth_ref(LOSidx)),radius*(ones(size(LOSidx,1),1)-elevation_ref(LOSidx)./90).*cosd(azimuth_ref(LOSidx)),'.','MarkerEdgeColor',[0.99 0.99 0.99]);
    plot(radius*(ones(size(NLOSidx,1),1)-elevation_ref(NLOSidx)./90).*sind(azimuth_ref(NLOSidx)),radius*(ones(size(NLOSidx,1),1)-elevation_ref(NLOSidx)./90).*cosd(azimuth_ref(NLOSidx)),'.','MarkerEdgeColor',[0.5 0.5 0.5]);

    idflg=1;
    figflg=1;
    flag=1;

    hlinex=[-0.9 0.9]; hliney=[0 0]; vlinex=[0 0]; vliney=[-0.9 0.9];
    arg=[0:100]*2*pi/100;
    % if figflg
    %    close
    % end
    plot(0.9*sin(arg),0.9*cos(arg),'k--',...
         0.8*sin(arg),0.8*cos(arg),'k--',...
         0.7*sin(arg),0.7*cos(arg),'k--',...
         0.6*sin(arg),0.6*cos(arg),'k--',...
         0.5*sin(arg),0.5*cos(arg),'k--',...
         0.4*sin(arg),0.4*cos(arg),'k--',...
         0.3*sin(arg),0.3*cos(arg),'k--',...
         0.2*sin(arg),0.2*cos(arg),'k--',...
         0.1*sin(arg),0.1*cos(arg),'k--',...
         hlinex,hliney,'k-',vlinex,vliney,'k-')
    axis('square')
    axis('off')
    % title(num2str(target_tow))
    % text(1.0,0,'E','FontSize',14)
    % text(-1.3,0,'W','FontSize',14)
    text(-0.1,0.85,'N','FontSize',14)
    % text(-0.1,-1,'S','FontSize',14);
    title(['Skymask - ',num2str(idt)]);

% --- plot nSV position ---
for idsvp = 1:1:size(SV_elaz{idt,1},1)
    if SV_elaz{idt,1}(idsvp,4) == 1
            plot(radius*(1-SV_elaz{idt,1}(idsvp,2)/90).*sind(SV_elaz{idt,1}(idsvp,3)),...
                 radius*(1-SV_elaz{idt,1}(idsvp,2)/90).*cosd(SV_elaz{idt,1}(idsvp,3)),...
                 '.','MarkerEdgeColor',[0,1,0],'MarkerSize',75);  
    elseif SV_elaz{idt,1}(idsvp,4) == 0
            plot(radius*(1-SV_elaz{idt,1}(idsvp,2)/90).*sind(SV_elaz{idt,1}(idsvp,3)),...
                 radius*(1-SV_elaz{idt,1}(idsvp,2)/90).*cosd(SV_elaz{idt,1}(idsvp,3)),...
                 '.','MarkerEdgeColor',[1,0,0],'MarkerSize',75);
    end
    if SV_elaz{idt,1}(idsvp,1)>=10
        text(radius*(1-SV_elaz{idt,1}(idsvp,2)/90).*sind(SV_elaz{idt,1}(idsvp,3))-0.08,...
             radius*(1-SV_elaz{idt,1}(idsvp,2)/90).*cosd(SV_elaz{idt,1}(idsvp,3)),...
             num2str(SV_elaz{idt,1}(idsvp,1)),'FontSize',15,'Color',[0,0,0],'FontWeight','Bold');
    else
        text(radius*(1-SV_elaz{idt,1}(idsvp,2)/90).*sind(SV_elaz{idt,1}(idsvp,3))-0.04,...
             radius*(1-SV_elaz{idt,1}(idsvp,2)/90).*cosd(SV_elaz{idt,1}(idsvp,3)),...
             num2str(SV_elaz{idt,1}(idsvp,1)),'FontSize',15,'Color',[0,0,0],'FontWeight','Bold');
    end
end
% saveas(h_gps,['result\skyplot\',num2str(idt)],'jpeg');
close(h_gps)
end
end

% check synchronization
% figure(2);
% hold on;
% plot(wls_llh(:,2),wls_llh(:,1),'k.-','MarkerSize',10);
% idt = 0;
% idt = idt + 1;
% plot(wls_llh(idt,2),wls_llh(idt,1),'r.','MarkerSize',25)
% plot(gt_lon_epoch(idt),gt_lat_epoch(idt),'b.','MarkerSize',25)
% plot_google_map('maptype', 'satellite');

disp(['2D Error -> ',num2str(mean(Error.WLS_2D)),'; 3D Error -> ',num2str(mean(Error.WLS_3D))]);












