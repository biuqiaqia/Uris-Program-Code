function [skymask, building, building_label] = generate_skymask_with_label(llh, building, accuracy)
% By NG Hoi-Fung, Ivan (ivannhf.ng@connect.polyu.hk)
% Input:
%     llh: position (lat, lon, alt)
%     kml: kml file path
%     accuracy: ~~ OPTIONAL ~~, the step to increase of building side, 0.005~0.015 should be fine
% Output:
%     skymask: skymask list (361x1 double)
    
    if ~exist('accuracy','var')
        accuracy = 0.015;
    end
    
    D2R = pi/180;
    R2D = 180/pi;
    
    skymask = zeros(361,1);
    
    xyz = llh2xyz([llh(1:2)*D2R, llh(3)]);
    
%     building = read_kml(kml);
    
    build_dist = zeros(length(building),1);
    for i = 1:length(building)
        if exist('distance') && exist('wgs84Ellipsoid') 
%             build_dist(i) = distance(llh(1),llh(2),building(i).center(1),building(i).center(2),wgs84Ellipsoid);
            build_dist(i) = m_idist(llh(2),llh(1),building(i).center(2),building(i).center(1),'wgs84');
%             m_idist(lon1,lat1,lon2,lat2,spheroid)
            continue;
        end
        if exist('build_dist')
            build_dist(i) = pdist([xyz;building(i).center_xyz]);
            continue;
        end
        build_dist(i) = sqrt(sum((xyz - building(i).center_xyz).^2));
    end
    [~,sort_idx] = sort(build_dist);
    building = building(sort_idx);
    
    building_label = ones(361,5)*-1; % 1:building index; 2:vertex index
    
    for i = 1:length(building)
        vllh = building(i).v_llh;
        if vllh(1,3) < llh(3); continue; end;
        vxyz = building(i).v_xyz;
        vaer = zeros(size(vllh,1),size(vllh,2));
        for j = 1:size(vllh,1)
            tllh = vllh(j,:);
            vaer(j,:) = llh2aer(tllh, llh);
        end
        [~,az_max_idx] = max(vaer(:,1));
        [~,az_min_idx] = min(vaer(:,1));

        if skymask(ceil(vaer(az_min_idx,1)):ceil(vaer(az_max_idx,1))) > max(vaer(:,2)); continue; end
        
        for j = 1:size(vllh,1)-1

            if skymask(ceil(min([vaer(j,1),vaer(j+1,1)])):ceil(max([vaer(j,1),vaer(j+1,1)]))) > max([vaer(j,2),vaer(j+1,2)]); continue; end

            xyz1 = vxyz(j,:);
            xyz2 = vxyz(j+1,:);
            vec = xyz2 - xyz1;
            unitv = normr(vec);
            mag = norm(vec);

            txyz = xyz1;

            while mag > norm(xyz1 - txyz)
                tllh = xyz2llh(txyz);
                tllh(1:2) = tllh(1:2)*R2D;
                taer = llh2aer(tllh, llh);

                if skymask(ceil(taer(1))) < taer(2)
                    skymask(ceil(taer(1))) = taer(2);
                    building_label(ceil(taer(1)), 1) = i;
                    building_label(ceil(taer(1)), 2) = j;
                    building_label(ceil(taer(1)), 3:5) = tllh;
                end
                txyz = txyz + unitv * accuracy;
            end
        end
    end
    if skymask(361) == 0; skymask(361) = skymask(1); end
end

function building = read_kml(kml)
    D2R = pi/180;
    R2D = 180/pi;
    
    fileID = fopen(kml,'r');
    line_count = 1;
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

    for i = 1:size(kml_str,1)
        C = strsplit(kml_str{i},{' ',',','>','<'});
        C = C(~cellfun(@isempty,C));
        C = str2double(C);
        C = C(~isnan(C));
        C = reshape(C,3,length(C)/3)';
        C = C(:,[2 1 3]);
        for j = 1:size(C,1)
            xyz(j,:) = llh2xyz([C(j,1:2)*D2R C(j,3)]);
        end
        building(i).v_llh = C;
        building(i).v_xyz = xyz;
        building(i).center = mean(C);
        building(i).center_xyz = mean(xyz);
        clear xyz, C;
    end
end

function pos_aer = llh2aer(pos_llh, org_llh)
% Input:
%   pos_llh: position llh
%   org_llh: origin llh
% Input:
%   pos_aer: position aer (radian)

D2R = pi/180;
R2D = 180/pi;

pos_xyz = llh2xyz([pos_llh(1:2)*D2R pos_llh(3)]);
org_xyz = llh2xyz([org_llh(1:2)*D2R org_llh(3)]);

pos_aer = xyz2aer(pos_xyz, org_xyz);
pos_aer(1:2) = pos_aer(1:2)*R2D;

end

function aer = xyz2aer( xyz, orgxyz )
%    INPUTS
%	xyz(1) = ECEF x-coordinate in meters
%	xyz(2) = ECEF y-coordinate in meters
%	xyz(3) = ECEF z-coordinate in meters
%
%	orgxyz(1) = ECEF x-coordinate of local origin in meters
%	orgxyz(2) = ECEF y-coordinate of local origin in meters
%	orgxyz(3) = ECEF z-coordinate of local origin in meters
%
%    OUTPUTS
%       aer:  Column vector
%		aer(1,1) = 'azimuth' - radian
%		aer(2,1) = 'elevation' - radian
%		aer(3,1) = 'radius' - meters

orgllh = xyz2llh(orgxyz);

[az,elev,slantRange] = ecef2aer(xyz(1), xyz(2), xyz(3), orgllh(1), orgllh(2), orgllh(3), wgs84Ellipsoid, 'radians');
aer = [az,elev,slantRange];

end

function xyz = llh2xyz(llh)
%LLH2XYZ  Convert from latitude, longitude and height
%         to ECEF cartesian coordinates.  WGS-84
%
%	xyz = LLH2XYZ(llh)	
%
%	llh(1) = latitude in radians
%	llh(2) = longitude in radians
%	llh(3) = height above ellipsoid in meters
%
%	xyz(1) = ECEF x-coordinate in meters
%	xyz(2) = ECEF y-coordinate in meters
%	xyz(3) = ECEF z-coordinate in meters

%	Reference: Understanding GPS: Principles and Applications,
%	           Elliott D. Kaplan, Editor, Artech House Publishers,
%	           Boston, 1996.
%
%	M. & S. Braasch 10-96
%	Copyright (c) 1996 by GPSoft
%	All Rights Reserved.

	phi = llh(1);
	lambda = llh(2);
	h = llh(3);

	a = 6378137.0000;	% earth semimajor axis in meters
	b = 6356752.3142;	% earth semiminor axis in meters	
	e = sqrt (1-(b/a).^2);

	sinphi = sin(phi);
	cosphi = cos(phi);
	coslam = cos(lambda);
	sinlam = sin(lambda);
	tan2phi = (tan(phi))^2;
	tmp = 1 - e*e;
	tmpden = sqrt( 1 + tmp*tan2phi );

	x = (a*coslam)/tmpden + h*coslam*cosphi;

	y = (a*sinlam)/tmpden + h*sinlam*cosphi;

	tmp2 = sqrt(1 - e*e*sinphi*sinphi);
	z = (a*tmp*sinphi)/tmp2 + h*sinphi;

	xyz(1) = x;
	xyz(2) = y;
	xyz(3) = z;
end

