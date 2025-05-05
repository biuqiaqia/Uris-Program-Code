% Instructions:
% 1. Please change the value of variables: kml_file, max_lat, min_lat, max_lon, min_lon
% 2. Program will construct the grid itself (also check the grid point inside building)
% 3. Define the start and end index you want to generate: st_i, ed_i
% 4. After finished the defined range on step 3, a temp file will be saved automatically on 'temp/<kml name>/<kml name>_<computer name>.mat'
% 5. Run 'skymask_merge_Ivan.m' to check the progress and save

clc;
clear;
close all;
warning('off','all');

D2R = pi/180;
R2D = 180/pi;

reso = 2;

% kml_file = 'kml/Frankfurt_P12.kml';
% max_lat = 50.1165302250929; 
% min_lat = 50.1131142149071;
% max_lon = 8.68132588316722;
% min_lon = 8.67601300483278;


% kml_file = 'kml/an_os.kml';
% % kml_file = 'gps_solution_TST/TST-Deep-Urban.kml';
% 
% % leftdown = [22.297331,114.172198;];
% % rightup = [22.300889,114.177524;];
% % max_lat = 22.305889; 
% % min_lat = 22.292031;
% % max_lon = 114.189524;
% % min_lon = 114.172198;
% 
% max_lat = 22.3032; 
% min_lat = 22.296943;
% max_lon = 114.1873;
% min_lon = 114.179042;
% idxname = 1;




% kml_file = 'kml/AN5.kml';
% max_lat = 22.296618; 
% min_lat = 22.294036;
% max_lon = 114.176610;
% min_lon = 114.1730501;

% kml_file = 'kml/AN6.kml';
% max_lat = 22.297339; 
% min_lat = 22.293899;
% max_lon = 114.171545;
% min_lon = 114.169460;
% kml_file = 'kml/AN8.kml';
% max_lat = 22.304062; 
% min_lat = 22.303031; % 22.3885780282722
% max_lon = 114.175972;
% min_lon = 114.173119; % 114.206513099381
% kml_file = 'kml/AN7.kml';
% max_lat = 22.303045; 
% min_lat = 22.300080;
% max_lon = 114.173182;
% min_lon = 114.170768;
% kml_file = 'kml/AN11.kml';
% max_lat = 22.297339; 
% min_lat = 22.293899;
% max_lon = 114.171545;
% min_lon = 114.169460;
% kml_file = 'kml/AN12.kml';
% max_lat = 22.306888; 
% min_lat = 22.305461;
% max_lon = 114.189345;
% min_lon = 114.187345;

% kml_file = 'kml/AN13.kml';
% max_lat = 22.306900; 
% min_lat = 22.305231;
% max_lon = 114.189497;
% min_lon = 114.187391;
% kml_file = 'kml/AN14.kml';
% max_lat = 22.300230; 
% min_lat = 22.295774;
% max_lon = 114.181304;
% min_lon = 114.176132;% kml_file = 'kml/AN9.kml';

% kml_file = 'kml/AN15.kml';
% max_lat = 22.305455; 
% min_lat = 22.302953;
% max_lon = 114.170766;
% min_lon = 114.169566;% kml_file = 'kml/AN9.kml';

% kml_file = 'kml/AN16.kml';
% max_lat = 22.300183; 
% min_lat = 22.296840;
% max_lon = 114.174481;
% min_lon = 114.170766;% kml_file = 'kml/AN9.kml';

% kml_file = 'kml/AN19.kml';
% max_lat = 22.305738; 
% min_lat = 22.302326;
% max_lon = 114.171200;
% min_lon = 114.169158;% kml_file = 'kml/AN9.kml';


% kml_file = 'kml/AN20.kml';
% % max_lat = 22.302602; 
% % min_lat = 22.296476;
% % max_lon = 114.180325;
% % min_lon = 114.175000;% kml_file = 'kml/AN9.kml';

% 
% kml_file = 'kml/AN21.kml';
% max_lat = 22.298013; 
% min_lat = 22.294056;
% max_lon = 114.176256;
% min_lon = 114.171422;% kml_file = 'kml/AN9.kml';


% kml_file = 'kml/AN22.kml';
% max_lat = 22.305772; 
% min_lat = 22.302541;
% max_lon = 114.170260;
% min_lon = 114.168733;% kml_file = 'kml/AN9.kml';

% kml_file = 'kml/AN23.kml';
% max_lat = 22.309055; 
% min_lat = 22.303300;
% max_lon = 114.172450;
% min_lon = 114.170900;

kml_file = 'kml/AN24.kml';
max_lat = 22.308000; 
min_lat = 22.305900;
max_lon = 114.187846;
min_lon = 114.184900;
% kml_file = 'kml/AN9.kml';
% max_lat = 22.303070;
% min_lat = 22.300012;
% max_lon = 114.177865;
% min_lon = 114.191992;

% kml_file = 'kml/AN10.kml';
% max_lat = 22.304203; 
% min_lat = 22.301203;
% max_lon = 114.191292;
% min_lon = 114.186549;
% kml_file = 'kml/Barcelona_P2_P3.kml';
% pos = [41.3824027800000,2.15425000000000];
% d = 200;
% lat2m = distance(pos(1),pos(2),pos(1)+1,pos(2),wgs84Ellipsoid); m2lat = 1/lat2m;
% lon2m = distance(pos(1),pos(2),pos(1),pos(2)+1,wgs84Ellipsoid); m2lon = 1/lon2m;
% max_lat = pos(1)+d*m2lat;
% min_lat = pos(1)-d*m2lat;
% max_lon = pos(2)+d*m2lon;
% min_lon = pos(2)-d*m2lon;

% kml_file = 'kml/Barcelona_P9_P10.kml';
% pos = [41.396205	2.153784];
% d = 220;
% lat2m = distance(pos(1),pos(2),pos(1)+1,pos(2),wgs84Ellipsoid); m2lat = 1/lat2m;
% lon2m = distance(pos(1),pos(2),pos(1),pos(2)+1,wgs84Ellipsoid); m2lon = 1/lon2m;
% max_lat = pos(1)+d*m2lat;
% min_lat = pos(1)-d*m2lat;
% max_lon = pos(2)+d*m2lon;
% min_lon = pos(2)-d*m2lon;

load('DTM_HK.mat');
fprintf('MSL data loaded.\n');

[kml_path, kml_name, kml_ext] = fileparts(kml_file);
temp_path = ['temp/',kml_name];

if ~exist(temp_path,'dir')
    mkdir(temp_path);
end

initial_file = [temp_path,'/build_',kml_name,'.mat'];
if ~exist(initial_file, 'file')
    fprintf('Reading %s ',kml_file);
    building = read_kml(kml_file);
    fprintf('... done\n');
    [lat_grid,lon_grid,alt_grid,skymask,height_grid,length_grid_x,length_grid_y] = construct_area(max_lat,max_lon,min_lat,min_lon,reso,building,DTM_HK,1);
%     [pointClould,pointClould_xyz] = construct_pointCloud(building, 0.004);
    [pointClould] = [];
    building_height = skymask;
    reflection_az = skymask;
    fprintf('x: %.2fm, y: %.2fm, total: %.2fm^2 (%d)/n', size(skymask,1)*reso, size(skymask,2)*reso, size(skymask,1)*reso*size(skymask,2)*reso, numel(skymask));
    save(initial_file,'lat_grid','lon_grid','alt_grid','skymask','building_height','reflection_az','building','max_lat','max_lon','min_lat','min_lon','reso','length_grid_x','length_grid_y','pointClould');
    fprintf('Initial file created.\n');
else
    load(initial_file);
    fprintf('Initial file loaded.\n');
end


temp_file = [temp_path,'/',kml_name,'_',getenv('computername')];
[~, temp_name, temp_ext] = fileparts(temp_file);
if exist([temp_file,'.mat'], 'file')
    load(temp_file);
    empty_idx = find(cellfun(@isempty, skymask));
    fprintf('%s mat file loaded, progress %d/%d, %.2f%%\n', temp_name, (numel(skymask)-length(empty_idx)), numel(skymask), (numel(skymask)-length(empty_idx))*100/numel(skymask));
else
    fprintf('Temp mat file not found\n');
end

S = cell(numel(skymask),1);
H = cell(numel(skymask),1);
A = cell(numel(skymask),1);
non_empty_idx = find(~cellfun(@isempty,skymask));
S(non_empty_idx) = skymask(non_empty_idx);
H(non_empty_idx) = building_height(non_empty_idx);
A(non_empty_idx) = reflection_az(non_empty_idx);

tic;
start_time = datetime('now');
po = gcp;
numWork = po.NumWorkers;

st_i = 1;
% st_i = [40000];
% ed_i = 89999;
ed_i = numel(skymask);
N = ed_i - st_i + 1;

D = parallel.pool.DataQueue;
h = waitbar(0, sprintf('Skymask generating, please wait... %d-%d',st_i,ed_i));
afterEach(D, @(x) nUpdateWaitbar(x,N,h));
global p;
p = 1;

for i = st_i:ed_i
    sTic = tic;
    
    if isempty(S{i})
        llh = [lat_grid(i),lon_grid(i),alt_grid(i)];
%         S{i} = generate_skymask(llh, kml_file, building, 0.01);
%         [S{i}, H{i}, A{i}] = generate_skymask_pointCloud(llh,pointClould);
        [S{i}, H{i}, A{i}] = generate_skymask(llh, kml_file, building, 1);
%         try
%             [S{i}, H{i}, A{i}] = generate_skymask_pointCloud(llh,pointClould);
%         catch
%             S{i} = generate_skymask(llh, kml_file, building, 0.01);
%         end
    end
    
    fprintf('Grid: %d/%d, %.2f, %.2f, %.2fs finished\n', i, numel(skymask), mean(S{i}), std(S{i}), toc/numWork);
    send(D, toc(sTic));
end
non_empty_idx = find(~cellfun(@isempty,S));
skymask(non_empty_idx) = S(non_empty_idx);
building_height(non_empty_idx) = H(non_empty_idx);
reflection_az(non_empty_idx) = A(non_empty_idx);

fprintf('%s: %d-%d Used time: %s\n',datestr(datetime('now'),'HH:MM:SS'),st_i,ed_i,datestr(abs(start_time - datetime('now')),'HH:MM:SS.FFF'));

% empty_idx = find(cellfun(@isempty, skymask));
% fprintf('Skymask progress %d/%d, %.2f%%\n', (numel(skymask)-length(empty_idx)), numel(skymask), (numel(skymask)-length(empty_idx))*100/numel(skymask));

if 0
    save([temp_file,'.mat'],'skymask','building_height','reflection_az');
    fprintf('Temp mat file saved.\n');
else
    save(['mat/',kml_name,'_0309.mat'],'skymask','length_grid_x','length_grid_y','lat_grid','lon_grid','alt_grid');
    fprintf('Whole mat file saved.\n');
end

fprintf('%s: Used time: %s\n',datestr(datetime('now'),'HH:MM:SS'),datestr(abs(start_time - datetime('now')),'HH:MM:SS.FFF'));
% disp([datestr(datetime('now'),'HH:MM:SS'), ': Used time: ', datestr(seconds(toc),'dd HH:MM:SS.FFF'), 's.']);


% function nUpdateWaitbar(x,N,h)
% global p
% waitbar(p/N, h);
% p = p + 1;
% fprintf('%.2f%% finished, %.2f s. \n',100*p/N, x);
% end
