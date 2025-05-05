function LOS_tag = check_LOS(endPoint_xyz,startPoint_xyz,bmodel,diff_pln_id)
% clc;
% clear;
% close all;

% load('data\mat\RT_diff_vis_debug.mat');
% % sv_xyz = diff_xyz;
% D2R = pi/180;
% R2D = 180/pi;
% sv_xyz = diff_xyz;
% bmodel = bmodel_lite0;

% check if a point is LOS in the model
LOS_tag = 1;
for idx = 1:size(bmodel,2)
    for jdx = 1:size(bmodel(idx).pln, 1)
        if idx == diff_pln_id(1) && ismember(jdx,diff_pln_id(:,2))
            continue;
        end
%        if idx == 214
%            fprintf('在第 %d 个数据处\n', idx);
%            keyboard;  % 暂停程序
%        end
        if LOS_tag == 1
            intersection = intersection_line_plane(startPoint_xyz, endPoint_xyz, bmodel(idx).pln(jdx,:));
%            fprintf('正在处理第 %d 个数据\n', idx);
            if check_lie_on_plane(bmodel,intersection,idx,jdx)==1 && (check_lie_between_2points(startPoint_xyz,endPoint_xyz,intersection)==1)
                LOS_tag = 0;
            end
        end
    end
end

end


% check
% D2R = pi/180;
% R2D = 180/pi;
% bmodel(idx).llh(jdx,:)
% int_llh = xyz2llh(intersection);
% int_llh = int_llh.*[R2D,R2D,1];
% norm(rcvr_xyz-intersection)


%%
% from Ivan.Ng
function [ intersection ] = intersection_line_plane( start_point, end_point, plane )
%   find the intersection between line and plane
%   input:
%       satellite_xyz (ECEF)
%       receiver_xyz (ECEF)
%       plane (plane equation constant: A, B, C, D)
line = end_point - start_point;
t = -(plane(4) + dot(start_point,plane(1:3))) / dot(line,plane(1:3));
intersection = start_point + t * line;
end

%%
% from Ivan.Ng
function [ lie_on_plane ] = check_lie_on_plane( buildings_struct, point, bulding_no, building_node )
D2R = pi/180;
R2D = 180/pi;
point1 = buildings_struct(bulding_no).xyz(building_node,:); % Vertex A
point1_llh = buildings_struct(bulding_no).llh(building_node,:);
point1_0lat = llh2xyz([point1_llh(1:2)*D2R 0.0]); % Vertex C
point2  = buildings_struct(bulding_no).xyz(building_node+1,:); % Vertex B
point2_llh  = buildings_struct(bulding_no).llh(building_node+1,:);
point2_0lat = llh2xyz([point2_llh(1:2)*D2R 0.0]); % Vertex D
vector_p1_p2 = point2 - point1; % vector AB
vector_p1_pt = point - point1; % vector AM
vector_p1_p1o = point1_0lat - point1; % vector AC
AMAB = dot(vector_p1_pt,vector_p1_p2);
ABAB = dot(vector_p1_p2,vector_p1_p2);
AMAC = dot(vector_p1_pt,vector_p1_p1o);
ACAC = dot(vector_p1_p1o,vector_p1_p1o);
lie_on_plane = logical((0 < AMAB) && (AMAB < ABAB) && (0 < AMAC) && (AMAC < ACAC));
end

%%
% from Ivan.Ng
function [ lies_on_flag ] = check_lie_between_2points( start_point, end_point, check_point )
%	Check the point lies on two points or not
%   output: 1-lies on, 0-not lies on 

lies_on_flag = 1; % default true;

vector_start_end = end_point - start_point;
dist_start_end = sqrt(sum((vector_start_end).^2));

vector_end_start = start_point - end_point;
dist_end_start = sqrt(sum((vector_end_start).^2));

vector_start_check = check_point - start_point;
dist_start_check = sqrt(sum((vector_start_check).^2));

vector_end_check = check_point - end_point;
dist_end_check = sqrt(sum((vector_end_check).^2));

lies_on_flag = logical((dot(vector_start_end, vector_start_check) >= 0) & (dist_start_end > dist_start_check));

end

%%












