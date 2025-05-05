% Input: kml file path
% Output: matlab structure for ray-tracing

function [building_struct] = kml2mat_lite(kml_file,el_mask,user_xyz)
[kml_str, building_str] = read_kml(kml_file);
[building_struct, building_count] = process_kml(building_str);

if ~isempty(el_mask)&&~isempty(user_xyz)
    btag = zeros(size(building_struct,2),1);
    for id_b = 1:1:size(building_struct,2)
        bxyz = mean(building_struct(id_b).xyz);
        bh = building_struct(id_b).llh(1,3);
        dist = norm(user_xyz-bxyz);
        bel = asind(bh/dist);
        if bel>el_mask
            btag(id_b,1)=1;
        end
    end
    building_struct = building_struct(btag==1);
end

end

function [ kml_str, building_str ] = read_kml( kml_file_str )
%   READ_KML file
%   Input: path of kml file
%   Output: string of kml file in each line

% kml_file_str = 'kml\TST_new_by_Ivan_20180820_org.kml';

% disp(['Reading ', kml_file_str, '....']);
kml_file = fopen(kml_file_str, 'r');
line_count = 1;
while (1)
    strln = fgetl(kml_file);
    
    kml_str{line_count,1} = strln;
    line_count = line_count + 1;
    
    if feof(kml_file)
        break;
    end
end
fclose(kml_file);
% disp(['Read completed. ']);

building_count = 1;
for idx = 1:length(kml_str)
    strln = kml_str{idx, 1};
    coor_str = strfind(strln, '<coordinates>');
    if ~isempty(coor_str)
        building_str{building_count, 1} = kml_str{idx+1, 1};
        building_count = building_count + 1;
    end
end


end

function [ buildings_struct, building_count ] = process_kml( building_str )
%   process kml string to buildings structure
%   input:
%           building_str: all building vertex in string(llh)
%   output:
%           buildings_struct:
%               building_llh,
%               building_xyz,
%               plane_equation(A,B,C,D)

D2R = pi/180;
R2D = 180/pi;

for building_count = 1:length(building_str)
    %     str = '114.176979622671,22.2980250151133,22 114.177047581559,22.2979667271138,22 114.176907503478,22.2978247684196,22 114.176848903067,22.2978750982748,22 114.176848025872,22.2979346444605,22 114.17693910858,22.298027107904,22 114.176979622671,22.2980250151133,22';
    str = building_str{building_count};
    coor = strsplit(str,{',', ' '},'CollapseDelimiters', false);
    coor = coor(~cellfun(@isempty,coor));
    %     coor = coor(1:end-1);
    no_coor = size(coor,2);
    building_cell = reshape(coor, 3, no_coor/3).';
    building_cell = building_cell(:,[2,1,3]);
    
    for idx = 1:length(building_cell)
        building_llh(idx,:) = [str2double(building_cell{idx,1}) str2double(building_cell{idx,2}) str2double(building_cell{idx,3})];
    end
    
    [lon, lat] = poly2cw(building_llh(:,2),building_llh(:,1));
    building_llh(:,1:2) = [lat lon];
    
    for idx = 1:length(building_cell)
        %         building_llh(idx,:) = [str2double(building_cell{idx,1}) str2double(building_cell{idx,2}) str2double(building_cell{idx,3})];
        building_xyz(idx,:) = llh2xyz([building_llh(idx,1:2)*D2R building_llh(idx,3)]);
        
        if(idx > 1) % find constant for plane
            pt1 = building_xyz(idx,:);
            pt2 = building_xyz(idx-1,:);
            pt3 = building_llh(idx-1,:);
            pt3(1:2) = pt3(1:2) * D2R;
            pt3(3) = 0.0;
            pt3 = llh2xyz(pt3);
            
            [a b c d] = find_plane_equation(pt1, pt2, pt3);
            building_equation(idx-1,:) = [a b c d];
        end
    end
    
    buildings_struct(building_count).node = length(building_cell);
    buildings_struct(building_count).llh = building_llh;
    buildings_struct(building_count).xyz = building_xyz;
    buildings_struct(building_count).plane = building_equation;
    
    clear building_llh;
    clear building_xyz;
    clear building_equation;
    
end

end

function [ a, b, c, d ] = find_plane_equation( point1, point2, point3 )
%FIND_PLANE_EQUATION Summary of this function goes here
%   Detailed explanation goes here

% vec_p1p2 = point2 - point1; % vector point 1 to point 2
% vec_p1p3 = point3 - point1; % vector point 1 to point 3
% 
% i = det([vec_p1p2(2:3);vec_p1p3(2:3)]);
% j = -det([vec_p1p2([1,3]);vec_p1p3([1,3])]);
% k = det([vec_p1p2(1:2);vec_p1p3(1:2)]);
% 
% a = i * point1(1);
% b = j * point1(2);
% c = k * point1(3);
% d = i * -1 * point1(1) + j * -1 * point1(2) + k * -1 * point1(3);

% 
% i = (vec_p1p2(2) * vec_p1p3(3)) - (vec_p1p2(3) * vec_p1p3(2));
% j = -((vec_p1p2(1) * vec_p1p3(3)) - (vec_p1p2(3) * vec_p1p3(1)));
% k = (vec_p1p2(1) * vec_p1p3(2)) - (vec_p1p2(2) * vec_p1p3(1));
% 
% a = i * point1(1);
% b = j * point1(2);
% c = k * point1(3);
% d = i * -1 * point1(1) + j * -1 * point1(2) + k * -1 * point1(3);

A = [point1;point2;point3];
A(:,1) = 1;
B = [point1;point2;point3];
B(:,2) = 1;
C = [point1;point2;point3];
C(:,3) = 1;
D = [point1;point2;point3];
D = D * -1;

a = det(A);
b = det(B);
c = det(C);
d = det(D);





end

