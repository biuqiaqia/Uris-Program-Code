clc;
clear;
close all;

D2R = pi/180;
R2D = 180/pi;

correct_skymask = 0; % Select correct skymask, shift 45deg

% building_model_path = 'kml\Frankfurt_P4.kml';
% building_model_path = 'kml\Frankfurt_P6.kml';
% building_model_path = 'kml\Frankfurt_P10.kml';
building_model_path = 'kml\Frankfurt_P12.kml';
% building_model_path = 'kml\Barcelona_P2_P3.kml';
% building_model_path = 'kml\Barcelona_P5.kml';
% building_model_path = 'kml\Barcelona_P9_P10.kml';

[filepath, filename, fileext] = fileparts(building_model_path);

db_path = ['temp\',filename];
% db_path = ['temp\'];
if ~exist(db_path, 'dir')
    disp('Directory not exist');
    return;
end

initial_file = [db_path,'\build_',filename,'.mat'];
load(initial_file);

basefolder = cd(db_path);
file_list = dir([filename,'*_*.mat']);

for i = 1:size(file_list,1)
    mat_file = ['mat',num2str(i)];
%     va = genvarname(mat_file);
    eval([mat_file, ' = load(file_list(',num2str(i),').name);']);
    eval([mat_file, '_idx = find(~cellfun(@isempty,',mat_file,'.skymask));']);
end

cd(basefolder);

fprintf('%d temp files found.\n', size(file_list,1));


for i = 1:size(file_list,1)
    eval(['skymask(mat',num2str(i),'_idx) = mat',num2str(i),'.skymask(mat',num2str(i),'_idx);']);
    eval(['building_height(mat',num2str(i),'_idx) = mat',num2str(i),'.building_height(mat',num2str(i),'_idx);']);
    eval(['reflection_az(mat',num2str(i),'_idx) = mat',num2str(i),'.reflection_az(mat',num2str(i),'_idx);']);
end
empty_idx = find(cellfun(@isempty,skymask));
empty_idx_h = find(cellfun(@isempty,building_height));
empty_idx_az = find(cellfun(@isempty,reflection_az));

%% Correct skymask
if correct_skymask
    for i = 1:numel(skymask)
        if isempty(skymask{i})|skymask{i}<=0; continue; end
        temp_sm = zeros(361,1);
        temp_fin_sm = zeros(361,1);
        temp_sm = skymask{i};
        temp_fin_sm(:) = temp_sm([90:-1:1,361:-1:91]);
        skymask(i) = {temp_fin_sm};
        clear temp_sm temp_fin_sm;
    end
end

%% Reshape table
building_idx = find(cellfun(@(x) (all(x<0 | x==90) & ~isempty(x)), skymask(:)));
complete_idx = find(cellfun(@(x) all(x>=0), skymask(:)));
incomplete_idx = find(cellfun(@(x) isempty(x), skymask(:)));
corner = [max(lat_grid(1,:)),max(lon_grid(:,1)); min(lat_grid(1,:)),min(lon_grid(:,1))];
fig = figure; hold on; axis equal;
plot(lon_grid(complete_idx),lat_grid(complete_idx),'o','MarkerEdgeColor','none','MarkerFaceColor','g');
plot(lon_grid(incomplete_idx),lat_grid(incomplete_idx),'o','MarkerEdgeColor','none','MarkerFaceColor','r');
plot(corner([3 3 4 4 3]),corner([2 1 1 2 2]),'b-');
plot(lon_grid(building_idx),lat_grid(building_idx),'o','MarkerEdgeColor','none','MarkerFaceColor',ones(1,3)*0.5);

fprintf('%.2f%% Skymask Fininshed: %d/%d\n', (numel(skymask)-length(empty_idx))*100/numel(skymask), numel(skymask)-length(empty_idx), numel(skymask));
if isempty(empty_idx) && isempty(empty_idx_h) && isempty(empty_idx_az)
    save(['mat\',filename,'_merged.mat'], 'skymask', 'building_height', 'reflection_az', 'lon_grid', 'lat_grid', 'alt_grid', 'length_grid_x', 'length_grid_y', 'reso', 'max_lat', 'max_lon', 'min_lat', 'min_lon');
    disp('Merged mat file saved.');
else
    disp('Mat not saved.');
end

