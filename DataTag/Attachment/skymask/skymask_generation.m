% ------------------------------------------------------------------------
% Main code for JION SDM
% State-of-art SDM based on Paul 2019 paper
%
% by GH.Zhang 2019/04/15
% guo-hao.zhang@connect.polyu.hk
% ------------------------------------------------------------------------

clc;
clear;
close all;

D2R = pi/180;
R2D = 180/pi;
m2lat = 1/110734;
m2lon = 1/103043;

% SUMO simulation
filename_obs = 'f9p_splitter_SDM_2';% Ublox
filename_nav = 'hksc0680';

idv = 6;

% Processing setup
data_type = 1; % 1-Rinex/2-mat
skymask_path = 'data/building_model/skymask_TST_ned.mat';
% load(skymask_path);
kml_file = 'gps_solution_TST/TST-Deep-Urban.kml';

leftdown = [22.297331,114.172198;];
rightup = [22.300889,114.176524;];
length_grid_x = floor(0.5*(rightup(:,1)-leftdown(:,1))/m2lat+0.5);
length_grid_y = floor(0.5*(rightup(:,2)-leftdown(:,2))/m2lon+0.5);

lat_grids = [];

for la = 1:length_grid_x
    lag = [];
    
    for lo = 1:length_grid_y
        lag(lo) = leftdown(:,1) + 2*la*m2lat;

    end;

    lat_grids(la,:) = lag;
end;
lat_grid = lat_grids;

lon_grids = [];
for lo = 1:length_grid_y
    log = [];
    
    for la = 1:length_grid_x
        log(la) = leftdown(:,2) + 2*lo*m2lon;

    end;

    lon_grids(:,lo) = log;
end;    
lon_grid = lon_grids;
    
[building_struct_all] = read_kml(kml_file);
% generate_skymask_with_label([lat_grid(1,1),lon_grid(1,1),5], building_struct_all,0.015)    


skymask = []

for la = 1:length_grid_x  
    
    sm = []
    
    for lo = 1:length_grid_y

%         sm(lo,:) = [generate_skymask_with_label([lat_grid(lo,la),lon_grid(lo,la),5], building_struct_all,0.015)];   
        skymask{la,lo} = [generate_skymask_with_label([lat_grid(la,lo),lon_grid(la,lo),5], building_struct_all,0.015)]; 
        
    end;
    
    
    
end;


save('0309.mat','lat_grid','lon_grid','length_grid_x','length_grid_y','skymask');
    
    