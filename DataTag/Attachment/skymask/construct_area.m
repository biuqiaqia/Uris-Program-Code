function [lat_grid,lon_grid,alt_grid,skymask,height_grid,length_grid_x,length_grid_y] = construct_area(max_lat,max_lon,min_lat,min_lon,reso,building,DTM_HK,plotting)
    D2R = pi/180;
    R2D = 180/pi;
    
    mean_lat = mean([max_lat min_lat]);
%     mean_lon = mean([max_lon min_lon]);
%     lat2m = abs(distance(mean_lat+1, mean_lon, mean_lat, mean_lon, wgs84Ellipsoid));
%     lon2m = abs(distance(mean_lat, mean_lon+1, mean_lat, mean_lon, wgs84Ellipsoid));
%     m2lat = 1/lat2m;
%     m2lon = 1/lon2m;
    m2lat = 1/110734;
    m2lon = 1/103043;
    lat_ = min_lat:reso*m2lat:max_lat;
    lon_ = min_lon:reso*m2lon:max_lon;
    
    [lat_grid,lon_grid] = meshgrid(lat_,lon_);
    
    alt_grid = zeros(size(lat_grid));
    height_grid = zeros(size(lat_grid));
    skymask = cell(size(lat_grid));
    
    lat_grid = round(lat_grid,8);
    lon_grid = round(lon_grid,8);
    alt_grid = round(alt_grid,4);
    
    length_grid_x = size(lat_grid,1);
    length_grid_y = size(lat_grid,2);
    
    bVertices = cell2mat(cellfun(@(x,y) [repmat(x,size(y,1),1),y],num2cell(1:length(building),1)',{building(:).v_llh;}','UniformOutput',false));
    
    for i = 1:numel(lat_grid)
        try
            alt_grid(i) = Get_HK_MSL(lat_grid(i),lon_grid(i),DTM_HK);
        catch
            alt_grid(i) = 20.0;
%             try
%                 alt_grid(i) = getElevations(lat_grid(i), lon_grid(i), 'key', 'AIzaSyAHE7Klqyxg2SPMExl-o41BXY9rM08jPks');
%             catch
%                 alt_grid(i) = 0.0;
%             end
        end
        vbuilding = unique(bVertices(distance(bVertices(:,2),bVertices(:,3),lat_grid(i),lon_grid(i),wgs84Ellipsoid)<150, 1));
        bFiltered = building(vbuilding);
        for j = 1:length(bFiltered)
            vllh = bFiltered(j).v_llh;
            in = inpolygon(lat_grid(i),lon_grid(i),vllh(:,1),vllh(:,2));
            if in == 1
                skymask(i) = {vllh(1,3)*-1};
                height_grid(i) = vllh(1,3)*-1;
                break;
            end
        end
    end
    
    try
        if exist('geoaxes') && exist('plotting','var') && plotting == 1
            plot_map(max_lat,max_lon,min_lat,min_lon,building);
        end
    catch
    end
end

function plot_map(max_lat,max_lon,min_lat,min_lon,building)
    area_user = [min_lat,min_lat,max_lat,max_lat,min_lat; min_lon,max_lon,max_lon,min_lon,min_lon;]';
    
    fig = figure(length(findobj('type','figure'))+1);
    gx = geoaxes;
    hold on;
    usr = geoplot(gca, area_user(:,1),area_user(:,2),'--', 'color', [0 0 0], 'linewidth', 1);
    if exist('building','var')
        center_list = [building.center];
        center_list = reshape(center_list,3,length(center_list)/3)';
        min_sim = min(center_list(:,3));
        max_sim = max(center_list(:,3));
        cmap = colormap((jet));
        cmap_max = size(cmap,1);
        caxis([min_sim, max_sim]);
        hcb = colorbar;
        title(hcb,'Building Height (m)');
        for i = 1:size(building,1)
            temp_llh_list = building(i).v_llh;
            height = temp_llh_list(1,3);
            cmap_idx = round((((cmap_max - 1) / (max_sim - min_sim)) * (height - min_sim)) + 1);
            geoplot(gca, temp_llh_list(:,1),temp_llh_list(:,2),'-', 'color', cmap(cmap_idx,:), 'linewidth', 2);
        end
    end
    legend([usr], {'Skymask area'});
    hold off;
    gx.Basemap = 'openstreetmap';
    % gx.Basemap = addCustomBasemap('openstreetmap','a.tile.openstreetmap.org');
    %     gx.MapCenter = [ground_truth_llh(1), ground_truth_llh(2)];
%     gx.MapCenter = [map_center(1), map_center(2)];
    % gx.ZoomLevel = 17;
    %     geolimits(gx,[map_center(1)-90*m2lat,map_center(1)+90*m2lat],[map_center(2)-90*m2lon,map_center(2)+90*m2lon]);
%     geolimits(gx,yrange,xrange);
    gx.Grid = 'off';
end



