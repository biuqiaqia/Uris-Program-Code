% Convert NMEA data to .mat data for ublox receiver
% 2018/01/26 by Guohao Zhang 

clc;
clear all;
close all;
format long g

exp_name = 'data\gnss_log_2022_11_04_14_53_01';
% filename_nav = 'data\eph\hksc008f.19n';
gt_waypoint = [22.297811,114.178289,4,1;...
               22.297558,114.177985,4,0];%LLH and epoch


%% Synchronize with rinex
obs_path = [exp_name,'.22o'];
rcvNmeaPathName = [exp_name,'.nmea'];
global weights 
weights = 0;
global FTABLE;
FTABLE=finv(0.9995,1,1:200)';
%enabled constellations
GPS_flag = 1;  GLO_flag = 1;  GAL_flag = 1;  BDS_flag = 1;  QZS_flag = 1; SBS_flag = 0;
GPS_col = 'b'; GLO_col = 'r'; GAL_col = 'g'; BDS_col = 'c'; QZS_col = 'm';

[constellations] = goGNSS.initConstellation(GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag);
constellations.nEnabledSat = constellations.nEnabledSat + 40;
nSatTot = constellations.nEnabledSat;

goGNSS;
%                 [Eph, iono] = load_RINEX_nav(filename_nav);            

                [pr1_R, ph1_R, pr2_R, ph2_R, dop1_R, dop2_R, snr1_R, snr2_R, ...
                 time_GPS, time_R, week_R, date_R, pos_R, interval, antoff_R, antmod_R] = ...
                 load_RINEX_obs(obs_path, constellations);
 
fprintf('RINEX files loading complete\n');
is_bias_tot=NaN(6,length(time_GPS));
[~, all_sow] = time2weektow(time_GPS);
disp(['---> Epoch Amount: ',num2str(size(all_sow,1))]);

%% time
GPSTimeS = all_sow(1);
GPSTimeE = all_sow(end);

%% Load rcv nmea and time table
nmeaGPSTIndexShift = 1;
nmeaLatIndexShift = 2; 
nmeaLonIndexShift = 4;
nmeaHeightMSIIndexShift = 9;
nmeaHeightGEOIndexShift = 11;
nmeaGPSTIndexShift
UTCGPSshift = 18;
secondsInaDay = 86400;
rcvNmea = textread(rcvNmeaPathName, '%s', 'whitespace', ',');
% rcvNmea = textread(rcvNmeaPathName, '%s', 'whitespace', ',','bufsize',4095*10);
rcvPosResult = zeros(1,3);
rcvPosResultCount = 0;
gpsTimeShift = fix(GPSTimeE/secondsInaDay)*secondsInaDay;
if GPSTimeE>=0 && GPSTimeS>=0 && GPSTimeE>GPSTimeS
    for rcvNmeaIdx = 1:size(rcvNmea,1)
        if ~isempty(rcvNmea{rcvNmeaIdx, 1})
        if size(rcvNmea{rcvNmeaIdx, 1},2)>1 && (rcvNmea{rcvNmeaIdx, 1}(1,1) == '$'||rcvNmea{rcvNmeaIdx, 1}(1,2) == '$')
            %rcvNmeaIdx
            time = str2double(rcvNmea{rcvNmeaIdx + nmeaGPSTIndexShift, 1});
            lat = str2double(rcvNmea{rcvNmeaIdx + nmeaLatIndexShift, 1});
            lon = str2double(rcvNmea{rcvNmeaIdx + nmeaLonIndexShift, 1});
            hour = fix(time/10000);
            minu = fix((time-hour*10000)/100);
            sec = fix(time-hour*10000-minu*100);
            UTCTime = hour*3600 + minu * 60 + sec
            GPSTime = UTCTime + UTCGPSshift + gpsTimeShift;       
            rcvNmea{34071,1}=0;
            height_msl = str2double(rcvNmea{rcvNmeaIdx + nmeaHeightMSIIndexShift, 1});
            height_geo = str2double(rcvNmea{rcvNmeaIdx + nmeaHeightGEOIndexShift, 1});
            
            if GPSTime>=GPSTimeS&&GPSTime<=GPSTimeE
                rcvPosResultCount = rcvPosResultCount+1;
                latDeg = fix(lat/100);
                latMin = (lat-latDeg*100);
                latInDegrees = latDeg+latMin/60.0;               
                
                lonDeg = fix(lon/100);
                lonMin = (lon-lonDeg*100);
                lonInDegrees = lonDeg+lonMin/60.0;
                
                rcvPosResult (rcvPosResultCount ,1) = GPSTime;
                rcvPosResult (rcvPosResultCount ,2) = latInDegrees;
                rcvPosResult (rcvPosResultCount ,3) = lonInDegrees;
                rcvPosResult (rcvPosResultCount ,4) = height_msl;
                rcvPosResult (rcvPosResultCount ,5) = height_geo;
                
            end
        end
        end        
    end
end

Id_select = [];
for ids = 1:1:size(rcvPosResult,1)
    if ~isnan(rcvPosResult(ids,2)) && rcvPosResult(ids,2)>20 && rcvPosResult(ids,2)<40 && ~isnan(rcvPosResult(ids,5))
        Id_select = [Id_select;ids];
    end
end

% idn = find(rcvPosResult(:,2) < 25);
rcvPosResult = rcvPosResult(Id_select,:);

% rcvPosResult(1,:)=[];
% dlmwrite([exp_name,'_nmea.csv'], rcvPosResult, 'precision', '%12.8f');

% ublox_data = csvread(['NMEA\',exp_name,'.csv']);

kmlwrite([exp_name,'_nmea.kml'],rcvPosResult(:,2),rcvPosResult(:,3),'Icon',...
    'http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png','IconScale',0.5,'Color',[1,0,0],'Name','  ');
% csvtokml('conventional\data_con.kml',ublox_data,'cyan_dot');

% rcvPosResult = [rcvPosResult(1:174,:);131574,rcvPosResult(174,2:5);rcvPosResult(175:end,:)];

return

if size(all_sow,1)==size(rcvPosResult,1)
    disp('Synchronization Succeed');
end

%% GT way point
% gt_waypoint
gt_lon_epoch = [];
gt_lat_epoch = [];
gt_waypoint(end,4)=size(all_sow,1);
for idwp = 1:1:size(gt_waypoint,1)-1
    gt_lon_epoch = [gt_lon_epoch,gt_waypoint(idwp,2):(gt_waypoint(idwp+1,2)-gt_waypoint(idwp,2))/(gt_waypoint(idwp+1,4)-gt_waypoint(idwp,4)):gt_waypoint(idwp+1,2)];
    gt_lat_epoch = [gt_lat_epoch,gt_waypoint(idwp,1):(gt_waypoint(idwp+1,1)-gt_waypoint(idwp,1))/(gt_waypoint(idwp+1,4)-gt_waypoint(idwp,4)):gt_waypoint(idwp+1,1)];
    gt_lon_epoch(end) = [];
    gt_lat_epoch(end) = [];
end
gt_lon_epoch = [gt_lon_epoch,gt_waypoint(end,2)];
gt_lat_epoch = [gt_lat_epoch,gt_waypoint(end,1)];

%% figure
figure(1)
hold on;
plot(rcvPosResult(:,3),rcvPosResult(:,2),'r.','MarkerSize',10)
xlabel('latitude');
ylabel('longitude');
plot_google_map('maptype', 'roadmap');

figure(2)
hold on;
plot(rcvPosResult(:,3),rcvPosResult(:,2),'k.-','MarkerSize',10)
plot(gt_lon_epoch(:),gt_lat_epoch(:),'m.-','MarkerSize',10)
idt = 0;
idt = idt + 1;
plot(rcvPosResult(idt,3),rcvPosResult(idt,2),'r.','MarkerSize',25)
plot(gt_lon_epoch(idt),gt_lat_epoch(idt),'b.','MarkerSize',25)
xlabel('latitude');
ylabel('longitude');
plot_google_map('maptype', 'roadmap');



