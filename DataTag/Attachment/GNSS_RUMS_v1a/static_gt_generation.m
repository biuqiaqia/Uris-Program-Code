% static location generation
% for GNSS RUMS opensky analysis to ION2021

clc;
clear;
close all;
D2R = pi/180;
R2D = 180/pi;

type = 2;

switch type
    case 1
%% Begin
load('data\RT_sim\11-SW-5A-all0.mat');
lat_grid_rt = lat_grid;
lon_grid_rt = lon_grid;
load('data\model\skymask_TST_ned.mat');
load('data\model\DTM_HK.mat','DTM_HK');
range_offset = 20;%meter

max_lat_sm = max(max(lat_grid));
min_lat_sm = min(min(lat_grid));
max_lon_sm = max(max(lon_grid));
min_lon_sm = min(min(lon_grid));
max_lat = max(max(lat_grid_rt));
min_lat = min(min(lat_grid_rt));
max_lon = max(max(lon_grid_rt));
min_lon = min(min(lon_grid_rt));

temp_xyz1 = llh2xyz([max_lat,max_lon,0].*[D2R,D2R,1]);
offset_xyz1 = enu2xyz([-range_offset,-range_offset,0],temp_xyz1);
offset_llh1 = xyz2llh(offset_xyz1).*[R2D,R2D,1];

temp_xyz2 = llh2xyz([min_lat_sm,min_lon_sm,0].*[D2R,D2R,1]);
offset_xyz2 = enu2xyz([+range_offset,+range_offset,0],temp_xyz2);
offset_llh2 = xyz2llh(offset_xyz2).*[R2D,R2D,1];

vmax_lat = offset_llh1(1);
vmin_lat = offset_llh2(1);
vmax_lon = offset_llh1(2);
vmin_lon = offset_llh2(2);

%% Random Sampling
N_sample = 0;
gt_data = cell(5,4);
N_open = 0;
while N_sample<50 || N_open<10
    pos_lat = offset_llh2(1) + rand*(offset_llh1(1)-offset_llh2(1));
    pos_lon = offset_llh2(2) + rand*(offset_llh1(2)-offset_llh2(2));    
    
    [~,idr_lat]= min(abs(lat_grid(1,:)-pos_lat));
	[~,idr_lon]= min(abs(lon_grid(:,1)-pos_lon));
    pos_skymask = skymask{idr_lon,idr_lat};
    
    if sum(pos_skymask) ~= 32490
        N_sample = N_sample + 1;
        H_MSL = Get_HK_MSL(pos_lat,pos_lon,DTM_HK);
        gt_data{N_sample,1} = [pos_lat,pos_lon,H_MSL,0,0,0];
        gt_data{N_sample,2} = 2;
        gt_data{N_sample,3} = pos_skymask;
        gt_data{N_sample,4} = mean(pos_skymask./90);
        if gt_data{N_sample,4}<=0.15
            N_open = N_open + 1;
        end
    end
end

%% Plot
figure(1)
hold on;
plot([min_lon,min_lon,max_lon,max_lon,min_lon],[min_lat,max_lat,max_lat,min_lat,min_lat],'k--','LineWidth',2);
plot([min_lon_sm,min_lon_sm,max_lon_sm,max_lon_sm,min_lon_sm],[min_lat_sm,max_lat_sm,max_lat_sm,min_lat_sm,min_lat_sm],'b--','LineWidth',2);
xlabel('Longitude (degree)');
ylabel('Latitude (degree)');
hBase = plot_openstreetmap('Alpha',1,'Scale',50000,'BaseUrl', "http://a.tile.openstreetmap.org");
plot([min_lon,min_lon,max_lon,max_lon,min_lon],[min_lat,max_lat,max_lat,min_lat,min_lat],'k--','LineWidth',2);
plot([min_lon_sm,min_lon_sm,max_lon_sm,max_lon_sm,min_lon_sm],[min_lat_sm,max_lat_sm,max_lat_sm,min_lat_sm,min_lat_sm],'b--','LineWidth',2);
plot([vmin_lon,vmin_lon,vmax_lon,vmax_lon,vmin_lon],[vmin_lat,vmax_lat,vmax_lat,vmin_lat,vmin_lat],'r--','LineWidth',2);

f_pos = figure(2);
hold on;
plot([vmin_lon,vmin_lon,vmax_lon,vmax_lon,vmin_lon],[vmin_lat,vmax_lat,vmax_lat,vmin_lat,vmin_lat],'r--','LineWidth',2);
xlabel('Longitude (degree)');
ylabel('Latitude (degree)');
hBase = plot_openstreetmap('Alpha',1,'Scale',50000,'BaseUrl', "http://a.tile.openstreetmap.org");
plot([vmin_lon,vmin_lon,vmax_lon,vmax_lon,vmin_lon],[vmin_lat,vmax_lat,vmax_lat,vmin_lat,vmin_lat],'r--','LineWidth',2);
for idv = 1:1:size(gt_data,1)
    if gt_data{idv,4}<=0.15
        plot(gt_data{idv,1}(1,2),gt_data{idv,1}(1,1),'ko','MarkerFaceColor','g','MarkerSize',10);
    elseif gt_data{idv,4}<=0.3
        plot(gt_data{idv,1}(1,2),gt_data{idv,1}(1,1),'ko','MarkerFaceColor','c','MarkerSize',10);
    elseif gt_data{idv,4}<=0.5
        plot(gt_data{idv,1}(1,2),gt_data{idv,1}(1,1),'ko','MarkerFaceColor','b','MarkerSize',10);
    else
        plot(gt_data{idv,1}(1,2),gt_data{idv,1}(1,1),'ko','MarkerFaceColor','r','MarkerSize',10);
    end
    text(gt_data{idv,1}(1,2),gt_data{idv,1}(1,1),num2str(idv),'FontSize',15,'Color',[0,0,0],'FontWeight','Bold');
end

figure(3);
hold on;
plot([vmin_lon,vmin_lon,vmax_lon,vmax_lon,vmin_lon],[vmin_lat,vmax_lat,vmax_lat,vmin_lat,vmin_lat],'r--','LineWidth',2);
xlabel('Longitude (degree)');
ylabel('Latitude (degree)');
hBase = plot_openstreetmap('Alpha',1,'Scale',50000,'BaseUrl', "http://a.tile.openstreetmap.org");
plot([vmin_lon,vmin_lon,vmax_lon,vmax_lon,vmin_lon],[vmin_lat,vmax_lat,vmax_lat,vmin_lat,vmin_lat],'r--','LineWidth',2);
for idv = 1:1:50%size(gt_data,1)
    if gt_data{idv,4}<=0.15
        plot(gt_data{idv,1}(1,2),gt_data{idv,1}(1,1),'ko','MarkerFaceColor','g','MarkerSize',10);
    elseif gt_data{idv,4}<=0.53
        plot(gt_data{idv,1}(1,2),gt_data{idv,1}(1,1),'ko','MarkerFaceColor','b','MarkerSize',10);
    else
        plot(gt_data{idv,1}(1,2),gt_data{idv,1}(1,1),'ko','MarkerFaceColor','r','MarkerSize',10);
    end
end

%% Saving position
duration = 180; % seconds
for idv = 1:1:size(gt_data,1)
    gt_data{idv,1} = ones(180,1).*gt_data{idv,1};
end

%% Manually add correlated pos
corr_pos_A = [22.298683,114.178784;...
              22.298652,114.178814;...
              22.298619,114.178846;...
              22.298586,114.178877;...
              22.298554,114.178908];
          
corr_pos_B = [22.299553,114.178281;...
              22.299515,114.178316;...
              22.299475,114.178355;...
              22.299434,114.178395;...
              22.299393,114.178436;...
              22.299358,114.178397;...
              22.299322,114.178356;...
              22.299286,114.178313;...
              22.299249,114.178271];

corr_pos = [corr_pos_A;corr_pos_B];
id_rand = size(gt_data,1);
for idc = 1:1:size(corr_pos,1)
    H_MSL = Get_HK_MSL(corr_pos(idc,1),corr_pos(idc,2),DTM_HK);
    gt_data{id_rand+idc,1} = ones(180,1).*[corr_pos(idc,1),corr_pos(idc,2),H_MSL,0,0,0];
    gt_data{id_rand+idc,2} = 2;
    
    [~,idr_lat]= min(abs(lat_grid(1,:)-corr_pos(idc,1)));
	[~,idr_lon]= min(abs(lon_grid(:,1)-corr_pos(idc,2)));
    pos_skymask = skymask{idr_lon,idr_lat};
    gt_data{id_rand+idc,3} = pos_skymask;
    gt_data{id_rand+idc,4} = mean(pos_skymask./90);
    
end
% save('data\gt_pos\static_test_2.mat','gt_data');%static pos file

%% New 10m spatial correlated agents
center_llh = [22.299108,114.178249,5];
center_xyz = llh2xyz([22.299108,114.178249,5].*[D2R,D2R,1]);
b1_xyz = llh2xyz([22.299600,114.178489,5].*[D2R,D2R,1]);
b2_xyz = llh2xyz([22.299886,114.178210,5].*[D2R,D2R,1]);
[b_az,~] = topocent(b2_xyz,b1_xyz);
dn = 10*cosd(b_az);
de = 10*sind(b_az);

set1 = [de*[-10:-1]',dn*[-10:-1]',zeros(10,1)];
set2 = [de*[1:9]',dn*[1:9]',zeros(9,1)];
set3 = [(-dn)*[1:5]',de*[1:5]',zeros(5,1)];

set_all = [set1;zeros(1,3);set2;set3];

set_llh = 0*set_all;
set_xyz = 0*set_all;
for id = 1:size(set_all,1)
    set_xyz(id,:) = enu2xyz(set_all(id,:),center_xyz);
    set_llh(id,:) = xyz2llh(set_xyz(id,:)).*[R2D,R2D,1];
end

kmlwrite('data\gt_pos\spatial_co.kml',set_llh(:,1),set_llh(:,2),'Icon',...
    'http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png','IconScale',0.5,'Color',[0,0,1],'Name','  ');

id_rand = 0;
for idc = 1:1:size(set_llh,1)
    gt_data{id_rand+idc,1} = ones(180,1).*[set_llh(idc,1),set_llh(idc,2),set_llh(idc,3),0,0,0];
    gt_data{id_rand+idc,2} = 2;
    
    [~,idr_lat]= min(abs(lat_grid(1,:)-set_llh(idc,1)));
	[~,idr_lon]= min(abs(lon_grid(:,1)-set_llh(idc,2)));
    pos_skymask = skymask{idr_lon,idr_lat};
    gt_data{id_rand+idc,3} = pos_skymask;
    gt_data{id_rand+idc,4} = mean(pos_skymask./90);
    
end

save('data\gt_pos\static_test_sc.mat','gt_data');%static pos file


    case 2
%% static generation direct pos
load('data\TST_4grids_testcases.mat');
load('data\model\DTM_HK.mat','DTM_HK');
total_epoch = 3600;
% load('data\model\skymask_TST_ned.mat');

for id = 1:1:size(test_cases.p1,1)
    p1_pos(id,1:2) = [lat_grid(test_cases.p1(id,1),test_cases.p1(id,2)),lon_grid(test_cases.p1(id,1),test_cases.p1(id,2))];
    p2_pos(id,1:2) = [lat_grid(test_cases.p2(id,1),test_cases.p2(id,2)),lon_grid(test_cases.p2(id,1),test_cases.p2(id,2))];
end
agent_pos = [p1_pos(1,:);p2_pos(1,:);p1_pos(2,:);p2_pos(2,:);p1_pos(3,:);p2_pos(3,:);];

for idc = 1:1:size(agent_pos,1)
    H_MSL = Get_HK_MSL(agent_pos(idc,1),agent_pos(idc,2),DTM_HK);
    gt_data{idc,1} = ones(total_epoch,1).*[agent_pos(idc,1),agent_pos(idc,2),H_MSL,0,0,0];
    gt_data{idc,2} = 2;
    
%     [~,idr_lat]= min(abs(lat_grid(1,:)-agent_pos(idc,1)));
% 	[~,idr_lon]= min(abs(lon_grid(:,1)-agent_pos(idc,2)));
%     pos_skymask = skymask{idr_lon,idr_lat};
%     gt_data{idc,3} = pos_skymask;
%     gt_data{id_rand+idc,4} = mean(pos_skymask./90);
    
end



figure(1)
hold on;
plot(agent_pos(1,2),agent_pos(1,1),'ko','MarkerFaceColor','g','MarkerSize',10);
plot(agent_pos(2,2),agent_pos(2,1),'ko','MarkerFaceColor','g','MarkerSize',10);
plot(agent_pos(3,2),agent_pos(3,1),'ko','MarkerFaceColor','r','MarkerSize',10);
plot(agent_pos(4,2),agent_pos(4,1),'ko','MarkerFaceColor','r','MarkerSize',10);
plot(agent_pos(5,2),agent_pos(5,1),'ko','MarkerFaceColor','b','MarkerSize',10);
plot(agent_pos(6,2),agent_pos(6,1),'ko','MarkerFaceColor','b','MarkerSize',10);
xlabel('Longitude (degree)');
ylabel('Latitude (degree)');
hBase = plot_openstreetmap('Alpha',1,'Scale',50000,'BaseUrl', "http://a.tile.openstreetmap.org");
plot(agent_pos(1,2),agent_pos(1,1),'ko','MarkerFaceColor','m','MarkerSize',10);
% plot(agent_pos(2,2),agent_pos(2,1),'ko','MarkerFaceColor','g','MarkerSize',10);
% plot(agent_pos(3,2),agent_pos(3,1),'ko','MarkerFaceColor','r','MarkerSize',10);
% plot(agent_pos(4,2),agent_pos(4,1),'ko','MarkerFaceColor','r','MarkerSize',10);
plot(agent_pos(5,2),agent_pos(5,1),'ko','MarkerFaceColor','b','MarkerSize',10);
plot(agent_pos(6,2),agent_pos(6,1),'ko','MarkerFaceColor','b','MarkerSize',10);


end




























































