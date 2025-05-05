% ------------------------------------------------------------------------
% GNSS Realistic Urban Measurement Simulatior [GNSS-RUMS]
% reformat based on Remote Sensing paper 
%
% support GPS/BDS
% single-user validation (labeled by SIM_Rinex.mat)
% with ION2019 data
%
% by GH.Zhang 2020/09/08
% guo-hao.zhang@connect.polyu.hk
% ------------------------------------------------------------------------

clc;
clear;
close all;
D2R = pi/180;
R2D = 180/pi;

%% testing data
load('data\gt_pos\gt_demo.mat','gt_data');% groundtruth data
eph_file = 'data\eph\hksc278a.21n';% ephemeris data
GPS_week = 2178;
sim_time = [172801,172801+size(gt_data{1,1},1)-1]; % starting gps_sec (change 172801 to your case)
building_kml = 'data\model\demo_HK_TST_fix.kml';% building model
tag_exc_BDSnew = 1;% exclude new launched BDS satellite to avoid bugs

%% simulation setup
goGNSS;
% ublox TCXO model
load('data\rcvr_para\ubx_m8t_para.mat');
% opensky-CNR model
load('data\rcvr_para\ubx_m8t_cnr_model_log_3.mat','fit_para');
% 3D building model elevation mask
bmodel_el_mask = 10;
% integer ambiguity
inam_range = 100;
N_am = round(rand(127,1)*inam_range);%fixed per SV
% recevier sampling frequency
rcvr_frq = 1;
% initial receiver clock bias
% dtR_ini = data{id_match,14}(1);%[by data or zero]
dtR_ini = 0;
% receiver clock drift for Doppler
dtRV_tag = 0;%[0-drift from model/1-initial drift from 1st epoch data/2-epochwise from data/3-from previous data 1st epoch]
% [VR_ini,dtRV_ini] = initial_rcvr_clock_drift_cal(data,id_match,gt_llh(1,:));
% load(ddt_file)
% cnr fitting type
para_cnr_fitting_type = 2;% 1-linear 2-logarithm fitting
% chip length
para_cl = goGNSS.V_LIGHT/1.023/1000000; 
% spacing between early&late correlator (chips) 
para_spacing = 1; 
% diffraction-simulation activation
diff_tag = 1; 
% diffraction occurence threshold (degree)
para_coeff_threshold = 0.1;
% diffraction occurence threshold (degree)
% para_diff_threshold = 10;
% reflection occurence threshold (degree)
% para_refl_threshold = 85;
% reflection occurence threshold (degree)
para_delay_threshold = 300;
% ephemeris loading 
[Eph, iono_para] = load_RINEX_nav(eph_file);
if tag_exc_BDSnew == 1
    Eph(:,Eph(30,:)>123)=[];
end
    
% GPS_time_s = weektow2time(GPS_week,sim_time(1),'G');
% GPS_time_e = weektow2time(GPS_week,sim_time(2),'G');

%% Simulation start
% idt = 15; %<====time control
for idn = 1:1:size(gt_data,1)
disp(['Simulation ID ==> ',num2str(idn),'/',num2str(size(gt_data,1))]);
GPS_time_s = weektow2time(GPS_week,sim_time(1)+gt_data{idn,1}(1,1)-1,'G');

gt_llh = gt_data{idn,1}(:,2:4);
end_ID = size(gt_llh,1);
ppm = ParforProgressbar(end_ID);
parfor idt = 1:end_ID %<== parallel here
% for idt = [467] %<== parallel here
% disp(['Simulation ID ==> ',num2str(idt)]);
timeR = GPS_time_s + (idt-1)*rcvr_frq;
timeR_ = mod(timeR,604800);
user_llh = gt_llh(idt,:);
user_pos = llh2xyz(user_llh.*[D2R,D2R,1]);%height = 4m(MSL)

%% Satellite Position Estimation
% rough estimation
eph_avail = Eph(30,:);
[XS, ~, ~, ~, ~, ~, ~] = satellite_positions(timeR, 20200*1e3*ones(1,length(eph_avail)), eph_avail, Eph, [], [], zeros(1,length(eph_avail)), zeros(1,length(eph_avail)), 0);
[az_all, el_all, ~] = topocent(user_pos, XS);
[~,ia,~] = unique(eph_avail); % negelct the duplicated ephs
idx_eph_vis_sat = find(el_all(ia)>0);
eph_avail_c = eph_avail(ia(idx_eph_vis_sat));
XS_avail = XS(ia(idx_eph_vis_sat),:);
pr_rough = zeros(size(XS_avail,1),1);
for ids = 1:1:size(XS_avail,1)
    pr_rough(ids,1) = norm(XS_avail(ids,:)-user_pos);
end

% refined estimation
pr_refined = zeros(size(XS_avail,1),1);
sat = eph_avail_c;
[XS, ~, ~, ~, ~, ~, ~] = satellite_positions(timeR, pr_rough, sat, Eph, [], [], zeros(1,length(sat)), zeros(1,length(sat)), 0);
for ids = 1:1:size(XS,1)
    pr_refined(ids,1) = norm(XS(ids,:)-user_pos);
end

% final estimation with dtR
delta_dt_r = (157.5 + sum(rcvr_para.err(idt,1:3)))/rcvr_para.f0;
dt_r = dtR_ini + delta_dt_r*(idt-1);

[XS, dtS, XS_tx, VS_tx, time_tx, no_eph, sys] = satellite_positions(timeR, pr_refined, sat, Eph, [], [], zeros(1,length(sat)), zeros(1,length(sat)), dt_r);
for ids = 1:1:size(XS,1)
    pr_refined(ids,1) = norm(XS(ids,:)-user_pos);
end

%SV velocity estimation (learned from RTKLIB)
tt = 1e-3;
[XS_tt, ~, ~, ~, ~, ~, ~] = satellite_positions(timeR+tt, pr_refined, sat, Eph, [], [], zeros(1,length(sat)), zeros(1,length(sat)), dt_r);        
VS = (XS_tt - XS) / tt;

%radio frequency
Eph_t = rt_find_eph(Eph, timeR, 127);
lambda = goGNSS.getGNSSWavelengths(Eph_t, 127);
lambda = lambda(sat);


%% General component simulation
% SV-user elevation/azimuth/true range
[az, el, range] = topocent(user_pos,XS);

% Receiver clock bias effect [unable to mimic the same value]
E_dtr = ones(size(sat,2),1).* (dt_r * goGNSS.V_LIGHT);

% Satellite clock bias
E_dts = dtS * goGNSS.V_LIGHT;

% Tropospheric delay
N = geoidheight(user_llh(1),user_llh(2));
e_h = user_llh(3) + N; %ellipsoid height
E_tropo = tropo_error_correction(el,e_h);

% Ionospheric delay
E_iono = iono_error_correction(user_llh(1),user_llh(2),az,el,timeR,iono_para,[]);

% Data-based receiver clock drift for Doppler
if dtRV_tag == 0
    dtRV_est_N = delta_dt_r;%GPS
    dtRV_est_B = delta_dt_r;%BDS
elseif dtRV_tag == 1
    dtRV_est_N = dtRV_ini(1);%GPS
    dtRV_est_B = dtRV_ini(1) + dtRV_ini(4);%BDS
elseif dtRV_tag == 2
    [~,dtRV_data] = initial_rcvr_clock_drift_cal(data,id_match-1+idt,user_llh);
    dtRV_est_N = dtRV_data(1);%GPS
    dtRV_est_B = dtRV_data(1) + dtRV_data(4);%BDS
end

% user velocity groundtruth
% if idt == 1
%     VR_gt = VR_ini';%initial from data
% %     VR_gt = [0,0,0];%initial from static
% else
%     user_pos_last = llh2xyz(gt_llh(idt-1,:).*[D2R,D2R,1]);
%     VR_gt = user_pos-user_pos_last;
% end

% user velocity groundtruth from SPAN/ENU
VR_gt = enu2xyz(gt_data{idn,1}(idt,5:7),user_pos)'-user_pos;
% user velocity groundtruth from static test
% VR_gt = [0,0,0];

%% Ray-tracing for GNSS signal paths
% building model reduction (old ray-tracing)
[bmodel_R_lite] = readkmlfile_lite(building_kml,bmodel_el_mask,user_pos);
[bmodel_D_lite] = kml2mat_lite(building_kml,bmodel_el_mask,user_pos);

% ray-tracing for reflection
[type, reflection, reflDist, reflALL, ~, ~, ~] = ray_tracing_refl_only(user_pos, XS, bmodel_R_lite, sat', timeR);

% ray-tracing for diffraction
diff_info = cell(size(XS,1),2);
for ids = 1:1:size(XS,1)
	% ray-tracing for diffraction
	diff_info{ids,1} = sat(ids);
	diff_info{ids,2} = ray_tracing_diff(XS(ids,:),user_pos,bmodel_D_lite,[]);
end

% ray-tracing result combined
RT_data = cell(size(XS,1),3);
for ids = 1:1:size(sat,2)
    RT_data{ids,1} = sat(ids);
    RT_SV_info = cell(1,6);
    id_ray = 1;
    % LOS signal
    if type(ids)==1 || type(ids)==2 
        RT_SV_info{id_ray,1} = 1;
        RT_SV_info{id_ray,2} = nan;
        RT_SV_info{id_ray,3} = 0;
        RT_SV_info{id_ray,4} = nan;
        RT_SV_info{id_ray,5} = nan;
        RT_SV_info{id_ray,6} = 1;
        id_ray = id_ray + 1;
    elseif type(ids)==4 || type(ids)==5 
        % bug
        disp('Bug exist...stop!');
%         return
    end
	% reflected signal
    for idr = 1:1:size(reflALL{ids,1},1)
        RT_SV_info{id_ray+idr-1,1} = 2;
        RT_SV_info{id_ray+idr-1,2} = reflALL{ids,1}(idr,1:3);
        RT_SV_info{id_ray+idr-1,3} = reflALL{ids,2}(idr,1);
        RT_SV_info{id_ray+idr-1,4} = nan;
        temp_inci_angle = incident_angle_cal(XS(ids,:),user_pos,reflALL{ids,1}(idr,1:3));
        [R_LHCP] = Reflection_coefficient_LHCP(reflALL{ids,1}(idr,1:3),XS(ids,:),user_pos,lambda(ids));
        RT_SV_info{id_ray+idr-1,5} = [0,temp_inci_angle];
        RT_SV_info{id_ray+idr-1,6} = abs(R_LHCP);
    end
    id_ray = id_ray + size(reflALL{ids,1},1);
    % diffracted signal   
    for idd = 1:1:size(diff_info{ids,2},1)
        RT_SV_info{id_ray+idd-1,1} = 3;
        RT_SV_info{id_ray+idd-1,2} = diff_info{ids,2}{idd,3};
        RT_SV_info{id_ray+idd-1,3} = diff_info{ids,2}{idd,4};
        RT_SV_info{id_ray+idd-1,4} = diff_info{ids,2}{idd,2};
        D_coeff = UTD_coefficient_sim(user_pos,XS(ids,:),diff_info{ids,2}{idd,3},diff_info{ids,2}{idd,2},lambda(ids));
        RT_SV_info{id_ray+idd-1,5} = [diff_info{ids,2}{idd,5}];
        RT_SV_info{id_ray+idd-1,6} = abs(D_coeff(1));
    end
    id_ray = id_ray + size(diff_info{ids,2},1);
    
    % ---exclude large diff and sort-by-delay
%     if id_ray>1
%         id_out = [];
%         for idv = 1:1:size(RT_SV_info,1)
%             if RT_SV_info{idv,5}(1)>0 && RT_SV_info{idv,5}(2)>=para_diff_threshold
%                 id_out = [id_out;idv];
%             end
% %             if RT_SV_info{idv,5}(1)==0 && RT_SV_info{idv,5}(2)>=para_refl_threshold
% %                 id_out = [id_out;idv];
% %             end
%         end
%         RT_data{ids,2} = sortrows(RT_SV_info(~ismember(1:size(RT_SV_info,1),id_out),:),3);
%     else
%         RT_data{ids,2} = [];
%     end
%     RT_data{ids,3} = RT_SV_info;
    
    % ---exclude large diff and sort-by-strength
%     if id_ray>1
%         id_out = [];
%         for idv = 1:1:size(RT_SV_info,1)
%             if RT_SV_info{idv,5}(1)>0 && RT_SV_info{idv,5}(2)>=para_diff_threshold
%                 id_out = [id_out;idv];
%             end
% %             if RT_SV_info{idv,5}(1)==0 && RT_SV_info{idv,5}(2)>=para_refl_threshold
% %                 id_out = [id_out;idv];
% %             end
%         end
%         RT_data{ids,2} = sortrows(RT_SV_info(~ismember(1:size(RT_SV_info,1),id_out),:),6,'descend');
%         if size(RT_data{ids,2},1)>1
%             RT_data{ids,2} = sortrows(RT_data{ids,2}(1:2,:),3);
%         end
%     else
%         RT_data{ids,2} = [];
%     end
%     RT_data{ids,3} = RT_SV_info;

    % ---over coefficient and sort by delay
    if id_ray>1
        id_out = [];
        for idv = 1:1:size(RT_SV_info,1)
            if RT_SV_info{idv,6}<para_coeff_threshold||RT_SV_info{idv,3}>para_delay_threshold
                id_out = [id_out;idv];
            end
        end
        RT_data{ids,2} = sortrows(RT_SV_info(~ismember(1:size(RT_SV_info,1),id_out),:),3);
    else
        RT_data{ids,2} = [];
    end
    RT_data{ids,3} = RT_SV_info;
end


%% signal type based measurement simulation
CNR_sim = nan(size(sat,2),6);
Epr_sim = zeros(size(sat,2),1);
pr_sim = zeros(size(sat,2),1);
phi_sim = zeros(size(sat,2),1);
Ephi_sim = zeros(size(sat,2),2);
delay_sim = nan(size(sat,2),6);
doppler_sim = zeros(size(sat,2),2);
R_coeff = zeros(size(sat,2),1);
sv_status = zeros(size(sat,2),2);

for ids = 1:1:size(sat,2) %Individual SV Simulation
    if isempty(RT_data{ids,2})
        % ===== NLOS without valid path =====
        CNR_sim(ids,1) = opensky_CN0_modeling(sat(ids),fit_para,el(ids),para_cnr_fitting_type);
        CNR_sim(ids,2:5) = nan;
        Ephi_sim(ids,1) = nan;
        delay_sim(ids,1:5) = [nan,nan,nan,nan,nan];
        if sat(ids)<=32
            ddt_user = dtRV_est_N;
        elseif sat(ids)>=86
            ddt_user = dtRV_est_B;
        end
        [doppler_sim(ids,1),doppler_sim(ids,2)] = Doppler_shift_simulation_lite(user_pos,VR_gt,XS(ids,:),VS(ids,:),ddt_user,lambda(ids),0,[]);
        % ===================================
    elseif size(RT_data{ids,2},1)==1 && RT_data{ids,2}{1,1}==1
        % =========== signle LOS ============
        sv_status(ids,1) = 1;
        CNR_sim(ids,1) = opensky_CN0_modeling(sat(ids),fit_para,el(ids),para_cnr_fitting_type);
        CNR_sim(ids,2:4) = nan;
        CNR_sim(ids,5) = CNR_sim(ids,1);
        Ephi_sim(ids,1) = 0;
        delay_sim(ids,1:5) = [0,nan,nan,nan,0];
        %Doppler
        if sat(ids)<=32
            ddt_user = dtRV_est_N;
        elseif sat(ids)>=86
            ddt_user = dtRV_est_B;
        end
        [doppler_sim(ids,1),doppler_sim(ids,2)] = Doppler_shift_simulation_lite(user_pos,VR_gt,XS(ids,:),VS(ids,:),ddt_user,lambda(ids),1,[]);
        % ===================================
    elseif size(RT_data{ids,2},1)==1 && RT_data{ids,2}{1,1}==2
        % =========== signle NLOS ===========
        sv_status(ids,1) = 2;
        [R_LHCP] = Reflection_coefficient_LHCP(RT_data{ids,2}{1,2},XS(ids,:),user_pos,lambda(ids));
        CNR0 = opensky_CN0_modeling(sat(ids),fit_para,el(ids),para_cnr_fitting_type);
        CNR_refl = Reflection_CNR_LHCP_lite(R_LHCP,CNR0,range(ids),RT_data{ids,2}{1,3});
        CNR_sim(ids,1) = CNR0;
        CNR_sim(ids,3:4) = nan;
        CNR_sim(ids,[2,5]) = CNR_refl;
        E_delay = RT_data{ids,2}{1,3};
        delay_sim(ids,1:5) = [nan,E_delay,nan,nan,E_delay];
        Ephi_sim(ids,1) = E_delay;
        R_coeff(ids,1) = R_LHCP;
        %Doppler
        if sat(ids)<=32
            ddt_user = dtRV_est_N;
        elseif sat(ids)>=86
            ddt_user = dtRV_est_B;
        end
        [doppler_sim(ids,1),doppler_sim(ids,2)] = Doppler_shift_simulation_lite(user_pos,VR_gt,XS(ids,:),VS(ids,:),ddt_user,lambda(ids),2,RT_data{ids,2}{1,2});
        % ===================================
    elseif size(RT_data{ids,2},1)==1 && RT_data{ids,2}{1,1}==3
        % =========== signle Diff ===========
        sv_status(ids,1) = 3;
        CNR0 = opensky_CN0_modeling(sat(ids),fit_para,el(ids),para_cnr_fitting_type);
        [CNR_diff,D_UTD] = UTD_GNSS_CNR_sim(user_pos,XS(ids,:),...
                                            RT_data{ids,2}{1,2},...%inter-point
                                            RT_data{ids,2}{1,4},...%o-face
                                            0,...%LOS_coeff
                                            0,...%LOS_delay
                                            lambda(ids),CNR0);
        CNR_sim(ids,1) = CNR0;
        CNR_sim(ids,[2,4]) = nan;
        CNR_sim(ids,[3,5]) = CNR_diff;
        E_delay = RT_data{ids,2}{1,3};
        Ephi_sim(ids,1) = E_delay;
        delay_sim(ids,1:5) = [nan,nan,E_delay,nan,E_delay];
        if sat(ids)<=32
            ddt_user = dtRV_est_N;
        elseif sat(ids)>=86
            ddt_user = dtRV_est_B;
        end
        [doppler_sim(ids,1),doppler_sim(ids,2)] = Doppler_shift_simulation_lite(user_pos,VR_gt,XS(ids,:),VS(ids,:),ddt_user,lambda(ids),2,RT_data{ids,2}{1,2});
        % ===================================
    elseif size(RT_data{ids,2},1)>1
        % ============ Multipath ============
        % ---C/N0 sim---
        CNR0 = opensky_CN0_modeling(sat(ids),fit_para,el(ids),para_cnr_fitting_type);
        Coeff_a = 0;
        Coeff_b = 0;
        % 1st signal
        if RT_data{ids,2}{1,1}==1
            %LOS
            Coeff_a = 1;
        elseif RT_data{ids,2}{1,1}==2
            %relfection
            Coeff_a = Reflection_coefficient_LHCP(RT_data{ids,2}{1,2},XS(ids,:),user_pos,lambda(ids));
        elseif RT_data{ids,2}{1,1}==3
            %diffraction
            D_coeff = UTD_coefficient_sim(user_pos,XS(ids,:),RT_data{ids,2}{1,2},RT_data{ids,2}{1,4},lambda(ids));
            Coeff_a = D_coeff(1);
        end
        type_a = RT_data{ids,2}{1,1};
        % 2nd signal
        if RT_data{ids,2}{2,1}==1
            %LOS
            Coeff_b = 1;
        elseif RT_data{ids,2}{2,1}==2
            %relfection
            Coeff_b = Reflection_coefficient_LHCP(RT_data{ids,2}{2,2},XS(ids,:),user_pos,lambda(ids));
        elseif RT_data{ids,2}{2,1}==3
            %diffraction
            D_coeff = UTD_coefficient_sim(user_pos,XS(ids,:),RT_data{ids,2}{2,2},RT_data{ids,2}{2,4},lambda(ids));
            Coeff_b = D_coeff(1);
        end
        type_b = RT_data{ids,2}{2,1};
        % superposition
        [CNR_MP,CNR_a,CNR_b] = CNR_superposition_cal(CNR0,lambda(ids),Coeff_a,RT_data{ids,2}{1,3},Coeff_b,RT_data{ids,2}{2,3});
        CNR_sim(ids,1) = CNR0;
        CNR_sim(ids,type_a) = CNR_a;
        if type_b~=type_a
            CNR_sim(ids,type_b) = CNR_b;
        else
            CNR_sim(ids,4) = CNR_b;%2nd CNR in 4th block
        end
        CNR_sim(ids,5) = CNR_MP;
        
        % ---pr sim---
        beta = -2*pi*(RT_data{ids,2}{2,3}-RT_data{ids,2}{1,3})/lambda(ids);%phase offset in rads
        [E_delay,Phi_delay] = multipath_pseudorange_delay(RT_data{ids,2}{1,3},RT_data{ids,2}{2,3},CNR_a,CNR_b,beta,para_spacing,para_cl);
        Ephi_sim(ids,1) = Phi_delay*lambda(ids)/pi;
        if type_b~=type_a
            delay_sim(ids,type_a) = RT_data{ids,2}{1,3};
            delay_sim(ids,type_b) = RT_data{ids,2}{2,3};
            delay_sim(ids,5) = E_delay;
        else
            delay_sim(ids,type_a) = RT_data{ids,2}{1,3};
            delay_sim(ids,4) = RT_data{ids,2}{2,3};%2nd Pr in 4th block
            delay_sim(ids,5) = E_delay;
        end
        
        % ---Doppler sim---
        if sat(ids)<=32
            ddt_user = dtRV_est_N;
        elseif sat(ids)>=86
            ddt_user = dtRV_est_B;
        end
        if type_a == 1
            % LOS available
            [doppler_sim(ids,1),doppler_sim(ids,2)] = Doppler_shift_simulation_lite(user_pos,VR_gt,XS(ids,:),VS(ids,:),ddt_user,lambda(ids),1,[]);
        elseif CNR_a>=CNR_b
            % signal a dominated
            [doppler_sim(ids,1),doppler_sim(ids,2)] = Doppler_shift_simulation_lite(user_pos,VR_gt,XS(ids,:),VS(ids,:),ddt_user,lambda(ids),2,RT_data{ids,2}{1,2});
        elseif CNR_a<CNR_b
            [doppler_sim(ids,1),doppler_sim(ids,2)] = Doppler_shift_simulation_lite(user_pos,VR_gt,XS(ids,:),VS(ids,:),ddt_user,lambda(ids),2,RT_data{ids,2}{2,2});
        end
        sv_status(ids,1) = type_a;
        sv_status(ids,2) = type_b;
        % ===================================
    end
end

%% General Noise Term
sigma_tropo = 0.12*1.001./sqrt(sind(el).^2+0.002001);
Mea_noise = zeros(size(XS,1),4);
for ids = 1:1:size(XS,1)
    if CNR_sim(ids,5)<5 % extreme weak signal neglected
        continue;
    end
    % Modeling Noise (SV_clock, Ephemeris, Iono, Tropo)
    sigma_model = sqrt(1.1^2 + 0.8^2 + 3.5^2 + sigma_tropo(ids)^2);
    % Tracking Loop Noise
    [sigma_DLL,sigma_PLL,sigma_FLL] = tracking_loop_error_modeling(CNR_sim(ids,5),lambda(ids),para_cl);
    Mea_noise(ids,:) = [randn*sigma_model,randn*sigma_DLL,randn*sigma_PLL,randn*sigma_FLL];
    % C/N0 noise added
    sigma_cnr = CNR_variance_modeling(CNR_sim(ids,5));
    CNR_sim(ids,6) = 10*log10(10^(CNR_sim(ids,5)/10)+randn*sigma_cnr);
    % Pr noise added
    delay_sim(ids,6) = delay_sim(ids,5) + Mea_noise(ids,1) + Mea_noise(ids,2);
    % Pr final sim
    pr_sim(ids,1) = range(ids) + delay_sim(ids,5) - E_dts(ids) + E_dtr(ids) + E_tropo(ids) + E_iono(ids) + Mea_noise(ids,1) + Mea_noise(ids,2);
    % Phi noise added
    Ephi_sim(ids,2) = Ephi_sim(ids,1) + Mea_noise(ids,3);
    % Phi final sim
    phi_sim(ids,1) = range(ids) + Ephi_sim(ids,1) - E_dts(ids) + E_dtr(ids) + E_tropo(ids) - E_iono(ids) + Mea_noise(ids,3) + lambda(ids)*N_am(sat(ids));
    % Doppler noise added
    doppler_sim(ids,3) = doppler_sim(ids,2) + Mea_noise(ids,4);
end


%% final result
RT_result{idt,idn} = RT_data;
CNR_result{idt,idn} = [sat',sv_status,CNR_sim];
PRE_result{idt,idn} = [sat',sv_status,el,az,delay_sim];
Phi_result{idt,idn} = [sat',sv_status,Ephi_sim,N_am(sat)];
Doppler_result{idt,idn} = [sat',sv_status,doppler_sim];
Mea_sim{idt,idn} = [sat',pr_sim,phi_sim,CNR_sim(:,6),doppler_sim(:,3)];
Noise_result{idt,idn} = Mea_noise;
PrCorrection{idt,idn} = - E_dts + E_tropo + E_iono;
SV_pos{idt,idn} = XS;
SV_vel{idt,idn} = [VS,lambda'];
ppm.increment();
end
delete(ppm);
end

%% Basic Simulation Output
Sim_data = cell(size(Mea_sim));
for idt = 1:size(Mea_sim,1) % for different epochs
    Sim_data{idt,1} = [GPS_week,sim_time(1)+idt-1];
    for idn = 1:size(Mea_sim,2) % for different agent
        if ~isempty(Mea_sim{idt,idn})
            Sim_data{idt,idn+1} = Mea_sim{idt,idn}(~isnan(Mea_sim{idt,idn}(:,2)),[1,2,4,5]);
        end
    end
end


%% Estimation Comparison
for idn = 1:size(gt_data,1)
    disp(['Positioning...',num2str(idn),'/',num2str(size(gt_data,1))]);
    for idt = 1:1:size(gt_data{idn,1},1)     
        % Least Squares Positioning
        temp_pos_gt = llh2xyz(gt_data{idn,1}(idt,2:4).*[D2R,D2R,1]);
        temp_SV_data_all = [Mea_sim{idt,idn}(:,1),Mea_sim{idt,idn}(:,2)-PrCorrection{idt,idn}(:,1),SV_pos{idt,idn},...
                            Doppler_result{idt,idn}(:,6),SV_vel{idt,idn},PRE_result{idt,idn}(:,4),CNR_result{idt,idn}(:,9)];

        temp_SV_data_all(isnan(temp_SV_data_all(:,2)),:)=[];
        [pos_sim_all(idt,:),~,~,~,W] = WLS_positioning(temp_SV_data_all(:,1:5),temp_SV_data_all(:,11),temp_SV_data_all(:,12),temp_pos_gt);
        llh_sim_all{idn,1}(idt,:) = xyz2llh(pos_sim_all(idt,:));
        llh_sim_all{idn,1}(idt,:) = llh_sim_all{idn,1}(idt,:).*[R2D,R2D,1];
        temp_pos_sim_all = llh2xyz([llh_sim_all{idn,1}(idt,1:2),gt_data{idn,1}(idt,4)].*[D2R,D2R,1]);
        Error_2D(idt,idn) = norm(temp_pos_gt-temp_pos_sim_all);
        denu_sim_all = xyz2enu(temp_pos_sim_all,temp_pos_gt);
        Error_E(idt,idn) = denu_sim_all(1);
        Error_N(idt,idn) = denu_sim_all(2);
        % Least Squares Velocity
        [VR,~,~,~,~] = LS_SA_code_Vel_Qmo_debug(pos_sim_all(idt,:)',...
                                                temp_SV_data_all(:,3:5),...
                                                temp_SV_data_all(:,7:9),...
                                                temp_SV_data_all(:,6),...
                                                [],temp_SV_data_all(:,11),...
                                                temp_SV_data_all(:,10),...
                                                ones(size(temp_SV_data_all,1),1)+(temp_SV_data_all(:,1)>=87).*3);
        V_enu{idn,1}(idt,:) = xyz2enu(temp_pos_gt+VR',temp_pos_gt);
        V_xyz{idn,1}(idt,:) = VR';
    end
    
    figure(idn)
    hold on;
    plot(gt_data{idn,1}(1,3),gt_data{idn,1}(1,2),'g.','MarkerSize',50);
    plot(gt_data{idn,1}(:,3),gt_data{idn,1}(:,2),'g-','LineWidth',3);
    plot(llh_sim_all{idn,1}(:,2),llh_sim_all{idn,1}(:,1),'ko','MarkerFaceColor','b','MarkerSize',7);
    xlabel('Longitude (degree)');
    ylabel('Latitude (degree)');
    hBase = plot_openstreetmap('Alpha',1,'Scale',50000,'BaseUrl', "http://a.tile.openstreetmap.org");
    plot(gt_data{idn,1}(1,3),gt_data{idn,1}(1,2),'g.','MarkerSize',50);
    plot(gt_data{idn,1}(:,3),gt_data{idn,1}(:,2),'g-','LineWidth',3);
    plot(llh_sim_all{idn,1}(:,2),llh_sim_all{idn,1}(:,1),'bo','MarkerFaceColor','b','MarkerSize',5);
end








































































