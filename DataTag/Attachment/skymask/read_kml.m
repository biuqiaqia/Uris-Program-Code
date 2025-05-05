function kml = readkml(kml_path)
    D2R = pi/180;
    R2D = 180/pi;
    
    % 读取KML文件
    fileID = fopen(kml_path, 'r');
    line_count = 1;
    kml_str = {};
    while ~feof(fileID)
        strln = fgetl(fileID);
        coor_str = strfind(strln, '<coordinates>');
        if ~isempty(coor_str)
            coor_str_ = strfind(strln, '</coordinates>');
            if ~isempty(coor_str_)
                kml_str{line_count,1} = strln;
            else
                strln = fgetl(fileID);
                kml_str{line_count,1} = strln;
            end
            line_count = line_count + 1;
        end
    end
    fclose(fileID);

    % 处理坐标信息
    for i = 1:size(kml_str,1)
        C = strsplit(kml_str{i},{' ',',','>','<'});
        C = C(~cellfun(@isempty,C));
        C = str2double(C);
        C = C(~isnan(C));
        C = reshape(C,3,length(C)/3)';
        C = C(:,[2 1 3]);
        
        % 第一种格式：包含经纬度、ECEF坐标和平面方程
        [x,y,z] = geodetic2ecef(wgs84Ellipsoid, C(:,1), C(:,2), zeros(size(C,1),1));
        e0 = [x,y,z];
        [x,y,z] = geodetic2ecef(wgs84Ellipsoid, C(:,1), C(:,2), C(:,3));
        e = [x,y,z];
        
        p = zeros(size(C,1)-1,4);
        for j = 1:size(C,1)-1
            pt1 = e(j+1,:);
            pt2 = e(j,:);
            pt3 = llh2xyz([C(j,1:2),0.0].*[D2R,D2R,1]);
            [a, b, c, d] = plane_eqn(pt1, pt2, pt3);
            p(j,:) = [a,b,c,d];
        end
        
        kml(i).llh = C;
        kml(i).xyz = e;
        kml(i).xyz0 = e0;
        kml(i).pln = p;
        
        % 第二种格式：包含经纬度、ECEF坐标和中心点
        %xyz = []; % 初始化 xyz 为空数组
        for j = 1:size(C,1)
            xyz(j,:) = llh2xyz([C(j,1:2)*D2R C(j,3)]);
        end
        kml(i).v_llh = C;
        kml(i).v_xyz = xyz;
        kml(i).center = mean(C);
        kml(i).center_xyz = mean(xyz);
        clear xyz;
    end
end

function [a, b, c, d] = plane_eqn(pt1, pt2, pt3)
    v1 = pt2 - pt1;
    v2 = pt3 - pt1;
    n = cross(v1, v2);
    a = n(1);
    b = n(2);
    c = n(3);
    d = -dot(n, pt1);
end

function xyz = llh2xyz(llh)
    lat = llh(1);
    lon = llh(2);
    alt = llh(3);
    a = 6378137.0; % WGS84 semi-major axis
    f = 1/298.257223563; % WGS84 flattening
    e2 = 2*f - f^2; % Square of eccentricity
    N = a / sqrt(1 - e2 * sin(lat)^2);
    x = (N + alt) * cos(lat) * cos(lon);
    y = (N + alt) * cos(lat) * sin(lon);
    z = (N * (1 - e2) + alt) * sin(lat);
    xyz = [x, y, z];
end