%code2
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

% Constants
D2R = pi/180;
R2D = 180/pi;

% Load building KML file
kml_file = 'file path.kml';
[building_struct_all] = read_kml(kml_file); % Ensure read_kml is implemented
% Read satellite data CSV file
obsname = '.obs file name';
satellite_filename = ['file path', obsname, '.csv']; % Path to satellite data file
satellite_data = readtable(satellite_filename); % Read satellite data
%%
% Extract necessary columns from CSV
azimuth = round(satellite_data.Azimuth); % Azimuth angle
elevation = satellite_data.Elevation; % Elevation angle
satellite_X = satellite_data.SaX; % Satellite X position (ECEF)
satellite_Y = satellite_data.SaY; % Satellite Y position (ECEF)
satellite_Z = satellite_data.SaZ; % Satellite Z position (ECEF)
XS = [satellite_X, satellite_Y, satellite_Z]; % Satellite ECEF coordinates
lambda = satellite_data.lambda(1);
%%
% Initialize visibility results
visibility = strings(size(azimuth)); % Store visibility results

% Define receiver position (ECEF coordinates)
receiver_lat = satellite_data.GT_Latitude(1); % Receiver latitude
receiver_lon = satellite_data.GT_Longitude(1); % Receiver longitude
receiver_height = satellite_data.GT_Height(1); % Receiver height
XR = llh2xyz([receiver_lat, receiver_lon, receiver_height] .* [D2R, D2R, 1]); % Receiver ECEF coordinates

% Call ray tracing function to analyze reflection signals
sat = satellite_data.PRN; % Satellite PRN
time_rx = satellite_data.GPS_Time; % GPS time

% Define grid points
lat_grid = receiver_lat;
lon_grid = receiver_lon;
h_grid = receiver_height;
GPS_Time = ['/'];

%%

% Generate skymask
skymask = [];
for idx = 1:length(lat_grid)  
    skymask{idx} = generate_skymask_with_label([lat_grid(idx), lon_grid(idx), h_grid(idx)], building_struct_all, 0.015); 
end

% Save skymask
save('MU_RL.mat', 'GPS_Time', 'skymask');

% Write skymask to file
fid = fopen('ACF_Loc1', 'w');
for idx = 1:length(lat_grid)  
    fprintf(fid, '%s', GPS_Time(idx));
    fprintf(fid, ','); 
    for j = 1:length(skymask{idx})
        fprintf(fid, '%f', skymask{idx}(j)); 
        fprintf(fid, ','); 
    end
    fprintf(fid, '\n'); 
end
fclose(fid);

% Read mask file ACF_Loc1
mask_filename = 'ACF_Loc1'; % Path to mask file
mask_data = fileread(mask_filename); % Read file content
mask_values = strsplit(mask_data, ','); % Split by comma
mask_elevation = str2double(mask_values(2:end)); % Extract elevation data
%}
%%



%%
%[Type, reflection, reflDist, diffraction, diffDist] = Recognize_type_REFLandDIFF(XR, XS, building_struct_all, sat, time_rx, azimuth);
[type, reflection, reflDist, diffraction, diffDist, reflALL, diffALL, D_COEFF, tooSmall, useTime] = ray_tracing_refl_and_diff(XR, XS, building_struct_all, sat, lambda, time_rx);

%%

% Determine visibility for each satellite
for i = 1:length(azimuth)
    az = azimuth(i); % Azimuth angle
    el = elevation(i); % Elevation angle
    % Check if azimuth is within valid range
    if az >= 1 && az <= length(mask_elevation)
        building_elevation = mask_elevation(az); % Building elevation at current azimuth
        
        % Determine visibility
        if el > building_elevation
            visibility(i) = 'Visible'; % Satellite is visible
        else
            visibility(i) = 'Blocked'; % Satellite is blocked
        end

    else
        fprintf('Satellite %d: Azimuth = %d° is invalid (must be between 1 and 360).\n', i, az);
    end
end
%}
%%


satellite_data.Visibility = visibility;
%satellite_data.ReflPoint = reflection;

satellite_data.ReflDist = reflDist;
%satellite_data.DifflPoint = diffraction;
satellite_data.DiffDist = diffDist;
satellite_data.D_Coeff = abs(D_COEFF);
satellite_data.Type = type;
satellite_data.tooSmall = tooSmall;

filename1 = ['DataTag/output/', obsname, '_tagDel.csv'];
NonDeletedDataType = satellite_data;
NonDeletedDataType = satellite_data(satellite_data.Type == 5, :);
NonDeletedDataType = satellite_data(satellite_data.Type == 6, :);
%NonDeletedDataType = satellite_data(satellite_data.tooSmall == 1, :);
writetable(NonDeletedDataType, filename1);

%satellite_data = satellite_data(satellite_data.tooSmall ~= 1, :);
satellite_data = satellite_data(satellite_data.Type ~= 5, :);
satellite_data = satellite_data(satellite_data.Type ~= 6, :);
satellite_data.tooSmall = [];

satellite_data = satellite_data(satellite_data.Elevation > 0, :);

% Save updated data to a new CSV file
output_filename = ['out put path', obsname, '_tag.csv']; % Output filename
writetable(satellite_data, output_filename);

% Print completion message
fprintf('Visibility and type analysis results saved to: %s\n', output_filename);
