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
% filename_obs = 'gnss_log_2023_08_18_16_10_14';% Ublox
% % load(['data\',filename_obs,'_nmea.mat']);
% % gt_waypoint = [rcvPosResult(1,2:4),1;...
% %                rcvPosResult(end,2:4),0];%LLH and epoch
% gt_waypoint = [22.304521586,114.17943231,3,1;...
%                22.3045215861,114.179432311,3.001,5000];
% filename_nav = 'data\hksc230i.23n';           
% ref_obs_path = 'data\hksc230i.23o';

% OPPO data
filename_obs = 'COM3___9600_231005_134614abc';% Ublox
% load(['data\',filename_obs,'_nmea.mat']);
% gt_waypoint = [rcvPosResult(1,2:4),1;...
%                rcvPosResult(end,2:4),0];%LLH and epoch
gt_waypoint = [22.304171,114.178915,12,1;...
               22.3041711,114.1789151,12.001,5000];
filename_nav = 'data\hksc278n.23n';           
ref_obs_path = 'data\hksc278n.23o';


% Processing Setup        
% file_path = 'data\building_model\demo_HK_TST_fix.kml';
% skymask_path = 'data\building_model\skymask_TST_ned';
ref_xyz = [-2414266.9197,5386768.9868,2407460.0314];%hksc
existing_skymask_tag = 1; % 0-no/ 1-exist/ 2-all open
skyplot_output = 0;
DD_pr_error_tag = 1;
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


%% Pseudorange Error from Double Difference
if DD_pr_error_tag == 1
    disp('-----calculating Pseudorange Error from DD-----');
    [Pr_error] = DD_pr_error_from_Ref(filename_nav,ref_obs_path,ref_xyz,data_all,gt_xyz);
    % check consistency
    if size(data_all,1) == size(Pr_error,1) 
        for idt = 1:size(Pr_error,1) 
            if size(data_all{idt,3},1) == 0
                check_diff(idt,1) = nan;
            elseif size(data_all{idt,3},1) == size(Pr_error{idt,1}(:,1),1)
                check_diff(idt,1) = 0;
            else
                check_diff(idt,1) = setdiff(data_all{idt,3},Pr_error{idt,1}(:,1));
            end
        end
    else
        disp('total epochs inconsistent');
    end


end

%% Plot
% WLS
ls_llh(ls_llh(:,1)==0,:)=nan;
wls_llh(wls_llh(:,1)==0,:)=nan;

figure(1);
hold on;
plot(ls_llh(:,2),ls_llh(:,1),'b.','MarkerSize',25);
plot(wls_llh(:,2),wls_llh(:,1),'r.','MarkerSize',25);
% plot(gt_epoch(:,2),gt_epoch(:,1),'g.-','MarkerSize',25);
hBase = plot_openstreetmap('Alpha',1,'Scale',50000,'BaseUrl', "http://a.tile.openstreetmap.org");
plot(ls_llh(:,2),ls_llh(:,1),'b.','MarkerSize',25);
plot(wls_llh(:,2),wls_llh(:,1),'r.','MarkerSize',25);
% plot(gt_epoch(:,2),gt_epoch(:,1),'g.-','MarkerSize',25);











