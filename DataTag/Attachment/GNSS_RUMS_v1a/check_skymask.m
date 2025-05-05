% skyplot & ray-tracing check
clc;
clear;
close all;
D2R = pi/180;
R2D = 180/pi;

%% Testing location & Building model
% Kowloon bay
user_llh = [22.319924,114.211178,5];
kml_file = 'data\model\Kowloon_bay.kml';
eph_file = 'data\eph\hksc297w.21n';
GPS_week = 2181;
sim_time = [81018,81027];

type = 2; % 1-skymask/2-ray-tracing

switch type
    case 1
%% generate skymask
        [building_struct_all] = read_kml(kml_file);
        [skymask,~,~] = generate_skymask_with_label(user_llh,building_struct_all);

        Fig = Skymask_plot_probability_NED([],skymask,[],1);
        
    case 2
%% generate ray-tracing
        %% simulation setup
        goGNSS;
        % 3D building model elevation mask
        bmodel_el_mask = 10;
        % recevier sampling frequency
        rcvr_frq = 1;
        % initial receiver clock bias
        dtR_ini = 0;
        % reflection occurence threshold (degree)
        para_delay_threshold = 300;
        % ephemeris loading 
        [Eph, iono_para] = load_RINEX_nav(eph_file);
        GPS_time_s = weektow2time(GPS_week,sim_time(1),'G');
        GPS_time_e = weektow2time(GPS_week,sim_time(2),'G');

        %% Simulation start
        idt = 10; %<====time control
        timeR = GPS_time_s + (idt-1)*rcvr_frq;
        timeR_ = mod(timeR,604800);
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
        dt_r = dtR_ini;
        [XS, dtS, XS_tx, VS_tx, time_tx, no_eph, sys] = satellite_positions(timeR, pr_refined, sat, Eph, [], [], zeros(1,length(sat)), zeros(1,length(sat)), dt_r);
        for ids = 1:1:size(XS,1)
            pr_refined(ids,1) = norm(XS(ids,:)-user_pos);
        end
        
        % ray-tracing
        [bmodel_R_lite] = readkmlfile_lite(kml_file,bmodel_el_mask,user_pos);
        [type, reflection, reflDist, reflALL, ~, ~, ~] = ray_tracing_refl_only(user_pos, XS, bmodel_R_lite, sat', timeR);
        
        RT_path_kml_out(sat',XS,user_pos,type,reflection,'result\RT_kml');
        
        [building_struct_all] = read_kml(kml_file);
        [skymask,~,~] = generate_skymask_with_label(user_llh,building_struct_all);
        
        % plot skymask [to be done]
        for idsv = 1 : length(XS)
            target_pos_xyz = [XS(idsv,1) XS(idsv,2) XS(idsv,3)];    
            target_pos_enu = xyz2enu(target_pos_xyz,user_pos);
            azimuth(idsv) = ((pi/2)-atan2(target_pos_enu(1),target_pos_enu(2)))*R2D;
            if (azimuth(idsv)<0)
                azimuth(idsv) = azimuth(idsv) + 360 ;
            end
            elevation(idsv) = ( R2D*atan(target_pos_enu(3)/norm(target_pos_enu(1:2))));
            if type(idsv) == 0 || type(idsv) == 3
                SV_vis(idsv) = 0;
            else
                SV_vis(idsv) = 1;
            end
        end
        azimuth_ned = azimuth.*(-1)+90;
        azimuth_ned(azimuth_ned<0) = azimuth_ned(azimuth_ned<0)+360;       
        SV_elaz = [sat',elevation',azimuth_ned',SV_vis'];

        Fig = Skymask_plot_probability_NED(SV_elaz,skymask,[],1);
end



















