function [type, reflection, reflDist, diffraction, diffDist, reflALL, diffALL, D_COEFF, tooSmall, useTime] = ray_tracing_refl_and_diff(XR, XS, kml, sat, lambda, time_rx)
% ========================================================================
% This version supports both reflection and diffraction calculations
% ========================================================================
% Input:
%   XR: receiver ECEF position (1 x 3)
%   XS: satellite ECEF position (n x 3)
%   kml: kml file path or processed structure of building model 
%   sat: satellite prn (OPTIONAL, if not correct/no input, will generate)
%   time_rx: time (OPTIONAL, will set 1 if no input)
% 
% Output:
%   type: signal type (-1: not valid; 1: LOS; 0: NLOS; 2: LOS with reflection; 3: NLOS with reflection; 4: NLOS with diffraction; 5: NLOS both refl & diff)
%   reflection: reflection ecef position with shortest delay
%   reflDist: reflection delay distance (meter)
%   diffraction: diffraction ecef position with shortest delay
%   diffDist: diffraction delay distance (meter)
%   reflALL: all possible reflection coordinate of each satellite (1st column: ecef; 2nd column: delay distance - meter)
%   diffALL: all possible diffraction coordinate of each satellite (1st column: ecef; 2nd column: delay distance - meter)
%   kml: processed structure of building model, can output and save for later use
%   data: formatted data
%       < time | prn | satellit lat | lon | alt | rover lat | lon | alt | type | reflection lat | lon | alt | diffraction lat | lon | alt >
%   useTime: time duration for ray-tracing (seconds)


rtTic = tic;
D2R = pi/180;
R2D = 180/pi;
warning('off','all');

% Load KML file if necessary
if ischar(kml)
    kml = readkml(kml);
end

% Handle empty inputs
if isempty(XR)
    type = [];
    reflection = [];
    reflDist = [];
    diffraction = [];
    diffDist = [];
    reflALL = [];
    diffALL = [];
    return; 
end

% Initialize variables
XRllh = xyz2llh(XR).*[R2D,R2D,1];
[az,el,r] = ecef2aer(XS(:,1),XS(:,2),XS(:,3),XRllh(1),XRllh(2),XRllh(3),wgs84Ellipsoid);


%%
%Debug



%%
% Find reflection points
[reflection, reflDist, reflALL, nlosIns] = find_reflection_points(XR, XS, kml, time_rx, sat);

% Find diffraction points
[isDiff, diffraction, diffDist, diffALL, D_COEFF, tooSmall] = find_diffraction_points(XS, XR, kml, lambda, time_rx, sat);

% Classify signal types
isRefl = reflection(:,1) ~= 0;
los = double(~logical(nlosIns));



%type = (~los&~isRefl&~isDiff).*0 + (los&~isRefl&~isDiff).*1 + (los&isRefl).*2 + (~los&isRefl&~isDiff).*3 + (~los&~isRefl&isDiff).*4 + (los&~isRefl&isDiff).*5 + (~los&isRefl&isDiff).*6;
%   .....NLOS..........................LOS.....................LOS_Refl...........NLOS_Refl........................LOS_Diff..................NLOS_Diff.............NLOS_Refl_Diff........
type = (~los&~isRefl&~isDiff).*0 + (los&~isRefl&~isDiff).*1 + (~los&isRefl&~isDiff).*2 + (~los&~isRefl&isDiff).*3 + (los&~isRefl&isDiff).*4 + (los&isRefl).*5 + (~los&isRefl&isDiff).*6;
%   .....NLOS..........................LOS.......................NLOS_Refl........................NLOS_Diff..................LOS_Diff.......


type(tooSmall == 1 & los == 0) = 0;   % 覆盖所有tooSmall=1的情况为NLOS
type(tooSmall == 1 & los == 1) = 1;   % 覆盖所有tooSmall=1的情况为LOS 


useTime = toc(rtTic); % get process time
end


























%% Sub-functions
function [reflection, reflDist, reflALL, nlosIns] = find_reflection_points(XR, XS, kml, time_rx, sat)
% Existing reflection point calculation logic from ray_tracing_refl_only
n = size(XS,1); % number of valid measurements
plnAll = vertcat(kml.pln);
plnN = num2cell([[1:length(kml)]',cell2mat(cellfun(@(x) size(x,1),[{kml.pln};]','UniformOutput',false))], 2);
plnNum = cell2mat(cellfun(@(x) [ones(x(2),1).*x(1), [1:1:x(2)]'], plnN, 'UniformOutput', false));
plnN = cell2mat(cellfun(@(x) x(2), plnN, 'UniformOutput', false));
xyzStart = cellfun(@(x) x(1:end-1,:),[{kml.xyz}]','UniformOutput',false); xyzStartAll = cell2mat(xyzStart);
xyzEnd = cellfun(@(x) x(2:end,:),[{kml.xyz}]','UniformOutput',false); xyzEndAll = cell2mat(xyzEnd);
xyz0Start = cellfun(@(x) x(1:end-1,:),[{kml.xyz0}]','UniformOutput',false); xyz0StartAll = cell2mat(xyz0Start);
xyz0End = cellfun(@(x) x(2:end,:),[{kml.xyz0}]','UniformOutput',false); xyz0EndAll = cell2mat(xyz0End);

nlosIns = zeros(size(XS,1),1);
ins = cell(size(XS,1),1);


for b = 1:length(kml)    
    mirror_pt = find_mirror_point(XR, kml(b).pln);
    roof = zeros(1,4);
    [roof(1), roof(2), roof(3), roof(4)] = plane_eqn(kml(b).xyz(1,:), kml(b).xyz(2,:), kml(b).xyz(3,:));
    for p = 1:size(mirror_pt,1)
        notBP = ~ismember(plnNum,[b,p],'rows');
        % find los/nlos
        ins_pt = intersection_line_plane(XR, XS, kml(b).pln(p,:));
        lie_on_plane = check_lie_on_plane(kml, ins_pt, b, p);  
        lies_on_flag = check_lie_between_2points(XR, XS, ins_pt);
        nlosIns = nlosIns + double(lie_on_plane&lies_on_flag);

        % find reflecting point
        ins_pt = intersection_line_plane(mirror_pt(p,:), XS, kml(b).pln(p,:));
        lie_on_plane = check_lie_on_plane(kml, ins_pt, b, p);
        idxRe = find(lie_on_plane);
        for r = 1:length(idxRe)
            sidx = idxRe(r);
            los_IS = ~check_line_plane_multiPln(ins_pt(sidx,:), XS(sidx,:), plnAll(notBP, :), xyzStartAll(notBP,:), xyz0StartAll(notBP,:), xyzEndAll(notBP,:)); % any intersection & satellite
            los_RI = ~check_line_plane_multiPln(XR, ins_pt(sidx,:), plnAll(notBP, :), xyzStartAll(notBP,:), xyz0StartAll(notBP,:), xyzEndAll(notBP,:)); % any intersection & satellite
            if all(los_IS&los_RI)
                ins{sidx} = [ins{sidx}; ins_pt(sidx,:)];
            end
        end
    end
end
idxRe = ~cell2mat(cellfun(@isempty, ins, 'UniformOutput', false));

%%
% calculate all delay distance
dist_RI = cellfun(@(x) sqrt(sum((XR-x).^2,2)), ins(idxRe), 'UniformOutput', false);
dist_IS = cellfun(@(x, y) sqrt(sum((x-y).^2,2)), ins(idxRe), num2cell(XS(idxRe,:),2), 'UniformOutput', false);
dist_RS = sqrt(sum((XR - XS(idxRe,:)).^2,2));
dist_refl = cellfun(@(x,y,z) x+y-z, dist_RI, dist_IS, num2cell(dist_RS,2), 'UniformOutput', false);
ins(idxRe,2) = dist_refl;

% exclude reflection delay with larger than 200m
inclIdx = cellfun(@(x) find(x < 200), ins(:,2), 'UniformOutput', false);
idxIncl = ~cell2mat(cellfun(@isempty, inclIdx, 'UniformOutput', false));
ins(idxIncl,1) = cellfun(@(x,y) x(y,:), ins(idxIncl,1), inclIdx(idxIncl), 'UniformOutput', false);
ins(idxIncl,2) = cellfun(@(x,y) x(y,:), ins(idxIncl,2), inclIdx(idxIncl), 'UniformOutput', false);

% find smallest delay as final reflection
[minDist, minIdx] = cellfun(@(x) min(x), ins(idxRe,2), 'UniformOutput', false);
reflection = zeros(n,3);
reflection(idxRe,:) = cell2mat(cellfun(@(x,y) x(y,:), ins(idxRe,1), minIdx, 'UniformOutput', false));
reflDist = zeros(n,1);
reflDist(idxRe,:) = cell2mat(cellfun(@(x,y) x(y,:), ins(idxRe,2), minIdx, 'UniformOutput', false));
reflALL = ins;
end

























function [isDiff, diffraction, diffDist, diffAll, D_COEFF, tooSmall] = find_diffraction_points(sv_xyz,rcvr_xyz,bmodel,lambda,time_rx, sat)
% Extracted diffraction point calculation logic from ray_tracing_diff_v3
D2R = pi/180;
R2D = 180/pi;

%initialization
diffraction = cell(size(time_rx));
diffDist = cell(size(time_rx));
diffAll = cell(size(time_rx, 1), 1);
isDiff = zeros(size(time_rx));
D_COEFF = zeros(size(time_rx));
tooSmall= zeros(size(time_rx));
for da = 1:1:size(sv_xyz, 1) % building-loop
    sv_xyz1 = sv_xyz(da, :);    % counter-clockwise model
    for idb = 1:1:size(bmodel,2)% building-loop
        center_xyz = bmodel(idb).center_xyz;
        clear bc_EAD;
        [bc_EAD(:,2),bc_EAD(:,1),bc_EAD(:,3)] = topocent(center_xyz,bmodel(idb).xyz);
        if sum((bc_EAD(1:size(bmodel(idb).xyz,1)-1,2)-bc_EAD(2:size(bmodel(idb).xyz,1),2))>0)<size(bmodel(idb).xyz,1)/2
            bmodel(idb).llh = flipud(bmodel(idb).llh);
            bmodel(idb).xyz = flipud(bmodel(idb).xyz);
            bmodel(idb).pln = flipud(bmodel(idb).pln);
        end
    end
    
    [sv_az,sv_el,sv_dist] = topocent(rcvr_xyz,sv_xyz1);
    
    diffall = cell(1,4);
    diff_edge_id = 1;
    for idb = 1:1:size(bmodel,2)% building-loop
        center_xyz = bmodel(idb).center_xyz;
        [center_AZ,~,~] = topocent(rcvr_xyz,center_xyz);
        if any(abs(center_AZ-sv_az)<90) || any(abs(center_AZ-sv_az)>270) % building in-front
            front_idc = [];
            for idc = 1:1:size(bmodel(idb).xyz,1)-1
                c1_xyz = bmodel(idb).xyz(idc,:);
                c2_xyz = bmodel(idb).xyz(idc+1,:);
                c3_xyz = (c1_xyz+c2_xyz)/2*0.99 + mean(bmodel(idb).xyz)*0.01;
                [c1_az,~,~] = topocent(rcvr_xyz,c1_xyz);
                [c2_az,~,~] = topocent(rcvr_xyz,c2_xyz);
    
                
                if (c1_az-c2_az)<0 || (c1_az-c2_az)>180 % front-edge
                    r1a = norm(c1_xyz-sv_xyz1);
                    r2a = norm(c2_xyz-sv_xyz1);
                    r1b = norm(c1_xyz-rcvr_xyz);
                    r2b = norm(c2_xyz-rcvr_xyz);
                    edge = norm(c1_xyz-c2_xyz);
    
                    % decent analysis
                    c1d_xyz = c1_xyz.*0.9999 + c2_xyz.*0.0001;
                    r_c1d = norm(c1d_xyz-sv_xyz1) + norm(c1d_xyz-rcvr_xyz);
                    c2d_xyz = c2_xyz.*0.9999 + c1_xyz.*0.0001;
                    r_c2d = norm(c2d_xyz-sv_xyz1) + norm(c2d_xyz-rcvr_xyz);
                    decent_tag = r_c1d<(r1a+r1b) & r_c2d<(r2a+r2b);
    
                    if decent_tag % diff-point is inside by convexity
                        diff_xyz = bi_section_nearest(c1_xyz,c2_xyz,rcvr_xyz,sv_xyz1);
                        % check blockage of incidence
                        diff_pln_id = [idb,idc];
                        LOS_tag1 = check_LOS(diff_xyz,rcvr_xyz,bmodel,diff_pln_id); % receiver-diffraction visible
                        LOS_tag2 = check_LOS(diff_xyz,sv_xyz1,bmodel,diff_pln_id); % diffraction-sateliite visible
                        if LOS_tag2 && LOS_tag1
                            diffall{diff_edge_id,1} = idb;% building ID
                            diffall{diff_edge_id,2} = [c1_xyz;c2_xyz;c3_xyz];% o-face points (RH rules)
                            corner_pos = [c1_xyz;c2_xyz;c3_xyz];% o-face points (RH rules)
                            diffall{diff_edge_id,3} = diff_xyz;
                            diffall{diff_edge_id,4} = norm(diff_xyz-sv_xyz1)+norm(diff_xyz-rcvr_xyz)-sv_dist;
                            [~,diff_el,~] = topocent(rcvr_xyz,diff_xyz);
                            diffall{diff_edge_id,5} = [1,abs(diff_el-sv_el)];
                            diff_edge_id = diff_edge_id+1;
                        end
                    end
                    front_idc = [front_idc;idc;idc+1];
                end
            end
            % Veritcal edge based on Horizontal edge
            for idv = 1:1:size(front_idc,1)
                % ----NEED LIMITATION----
                idc = front_idc(idv);
                if idc == 1
                    cL_xyz = bmodel(idb).xyz(end-1,:);
                    cR_xyz = bmodel(idb).xyz(idc+1,:);
                        
                elseif idc == size(bmodel(idb).xyz,1)
                    cL_xyz = bmodel(idb).xyz(idc-1,:);
                    cR_xyz = bmodel(idb).xyz(2,:);
                else
                    cL_xyz = bmodel(idb).xyz(idc-1,:);
                    cR_xyz = bmodel(idb).xyz(idc+1,:);
                end
                
                [cL_az,~,~] = topocent(bmodel(idb).xyz(idc,:),cL_xyz);
                [cR_az,~,~] = topocent(bmodel(idb).xyz(idc,:),cR_xyz);
                
                valid_ve = 0;
                if cL_az>cR_az
                    valid_ve = any(sv_az>cR_az) && any(sv_az<cL_az);
                else
                    valid_ve = any(sv_az<cL_az) || any(sv_az>cR_az);
                end
                if valid_ve~=1
                    continue
                end
                
                if rem(idv,2)
                    % left vertical edge
                    c1_xyz = llh2xyz(bmodel(idb).llh(idc,:).*[D2R,D2R,0]);%ground
                    c2_xyz = bmodel(idb).xyz(idc,:);%roof 
                    if idc == 1
                        c3_xyz = llh2xyz(bmodel(idb).llh(end-1,:).*[D2R,D2R,0]);%ground
                    else
                        c3_xyz = llh2xyz(bmodel(idb).llh(idc-1,:).*[D2R,D2R,0]);%ground
                    end
                else
                    % right vertical edge
                    c1_xyz = bmodel(idb).xyz(idc,:);%roof
                    c2_xyz = llh2xyz(bmodel(idb).llh(idc,:).*[D2R,D2R,0]);%ground
                    if idc == size(bmodel(idb).xyz,1)
                        c3_xyz = bmodel(idb).xyz(2,:);%roof
                    else
                        c3_xyz = bmodel(idb).xyz(idc+1,:);%roof
                    end
                end
    
                [edge_az,~]=topocent(rcvr_xyz,c1_xyz);
                
                % check angle inbetween requirement
                r1a = norm(c1_xyz-sv_xyz1);
                r2a = norm(c2_xyz-sv_xyz1);
                r1b = norm(c1_xyz-rcvr_xyz);
                r2b = norm(c2_xyz-rcvr_xyz);
                edge = norm(c1_xyz-c2_xyz);
                % corner 1 & 2
                beta1 = acosd((r1b^2 + edge^2 - r2b^2)/(2*r1b*edge));
                beta1_p = acosd((r1a^2 + edge^2 - r2a^2)/(2*r1a*edge));
                beta2 = acosd((r2b^2 + edge^2 - r1b^2)/(2*r2b*edge));
                beta2_p = acosd((r2a^2 + edge^2 - r1a^2)/(2*r2a*edge));                
                % supplementary angle
                beta1 = beta1 + 2*(beta1>90)*(90-beta1);
                beta1_p = beta1_p + 2*(beta1_p>90)*(90-beta1_p);
                beta2 = beta2 + 2*(beta2>90)*(90-beta2);
                beta2_p = beta2_p + 2*(beta2_p>90)*(90-beta2_p);
                % diff-point inside edge
                if (beta1>beta1_p)~=(beta2>beta2_p)
                    diff_xyz = bi_section_nearest(c1_xyz,c2_xyz,rcvr_xyz,sv_xyz1);
                    diff_dist = norm(diff_xyz-sv_xyz1)+norm(diff_xyz-rcvr_xyz)-sv_dist;
                    if idc == 1
                        diff_pln_id = [idb,idc;idb,size(bmodel(idb).xyz,1)-1];
                    elseif idc == size(bmodel(idb).xyz,1)
                        diff_pln_id = [idb,1;idb,idc-1];
                    else
                        diff_pln_id = [idb,idc-1;idb,idc];
                    end
    
                    % check blockage of incidence
                    LOS_tag1 = check_LOS(rcvr_xyz,diff_xyz,bmodel,diff_pln_id); % receiver-diffraction visible
                    LOS_tag2 = check_LOS(diff_xyz,sv_xyz1,bmodel,diff_pln_id); % diffraction-sateliite visible
                    if LOS_tag2 && LOS_tag1 && ~ismember(diff_dist,[diffall{:,4}])
                        diffall{diff_edge_id,1} = idb;% building ID
                        diffall{diff_edge_id,2} = [c1_xyz;c2_xyz;c3_xyz];% o-face points (RH rules)
                        corner_pos = [c1_xyz;c2_xyz;c3_xyz];% o-face points (RH rules)
                        diffall{diff_edge_id,3} = diff_xyz;
                        diffall{diff_edge_id,4} = diff_dist;
                        [diff_az,~,~] = topocent(rcvr_xyz,diff_xyz);
                        ddiff_az = abs(diff_az-sv_az);
                        if ddiff_az>180
                            ddiff_az = 360-ddiff_az;
                        end
                        diffall{diff_edge_id,5} = [2,ddiff_az];
                        diff_edge_id = diff_edge_id+1;
                    end
                end
                % ---end each vertical edge
            end
        end
    end
    
    disp(size(diffall)); 

    diffAll{da} = diffall;
    if all(cellfun(@isempty, diffall))
        diffDist{da} = 0;
        diffraction{da} = zeros(1,3);
        isDiff(da) = 0;
    else
        fprintf('diffpoint: %d', diff_xyz);
        D_coeff = UTD_coefficient_sim(rcvr_xyz,sv_xyz1,diff_xyz,corner_pos,lambda);
        D_COEFF(da) = D_coeff(1);
        if abs(D_COEFF(da)) > 0.1
            tooSmall(da) = 0;
        else
            tooSmall(da) = 1;
        end
        diffDist{da} = min(cell2mat(diffall(:, 4))); %存储路径差值
        diffraction{da} = diffall(:, 3);
        isDiff(da) = 1;
    end
    
end

end




















%% Helper functions
function mirror_pt = find_mirror_point(xyz, pln)
t = -(pln(:,4) + sum(xyz.*pln(:,1:3),2)) ./ (sum(pln(:,1:3).^2,2));
mirror_pt = xyz + 2.*t.*pln(:,1:3);
return;
end

function intersection = intersection_line_plane(start_point, end_point, plane)
line = end_point - start_point;
t = -(plane(4) + dot(start_point,plane(1:3),2)) ./ dot(line,repmat(plane(1:3),size(line,1),1),2);
intersection = start_point + t .* line;
end

function lie_on_plane = check_lie_on_plane(buildings_struct, point, bulding_no, building_node)
D2R = pi/180;
R2D = 180/pi;
n = size(point,1);
point1 = buildings_struct(bulding_no).xyz(building_node,:); % Vertex A
point1_llh = buildings_struct(bulding_no).llh(building_node,:);
point1_0lat = llh2xyz([point1_llh(1:2)*D2R 0.0]); % Vertex C
point2  = buildings_struct(bulding_no).xyz(building_node+1,:); % Vertex B
point2_llh  = buildings_struct(bulding_no).llh(building_node+1,:);
point2_0lat = llh2xyz([point2_llh(1:2)*D2R 0.0]); % Vertex D
vector_p1_p2 = point2 - point1; % vector AB
vector_p1_pt = point - point1; % vector AM
vector_p1_p1o = point1_0lat - point1; % vector AC
AMAB = dot(vector_p1_pt,repmat(vector_p1_p2,n,1),2);
ABAB = dot(vector_p1_p2,vector_p1_p2);
AMAC = dot(vector_p1_pt,repmat(vector_p1_p1o,n,1),2);
ACAC = dot(vector_p1_p1o,vector_p1_p1o);
lie_on_plane = logical((0 < AMAB) & (AMAB < ABAB) & (0 < AMAC) & (AMAC < ACAC));
end

function lies_on_flag = check_lie_between_2points(start_point, end_point, check_point)
lies_on_flag = 1; % default true;
vector_start_end = end_point - start_point;
dist_start_end = sqrt(sum((vector_start_end).^2, 2));
vector_end_start = start_point - end_point;
dist_end_start = sqrt(sum((vector_end_start).^2, 2));
vector_start_check = check_point - start_point;
dist_start_check = sqrt(sum((vector_start_check).^2, 2));
vector_end_check = check_point - end_point;
dist_end_check = sqrt(sum((vector_end_check).^2, 2));
try
    lies_on_flag = logical((dot(repmat(vector_start_end,size(check_point,1),1), vector_start_check, 2) >= 0) & (dist_start_end > dist_start_check)) ...
        & logical((dot(repmat(vector_end_start,size(check_point,1),1), vector_end_check, 2) >= 0) & (dist_end_start > dist_end_check));
catch
    lies_on_flag = logical((dot(vector_start_end, vector_start_check, 2) >= 0) & (dist_start_end > dist_start_check)) ...
        & logical((dot(vector_end_start, vector_end_check, 2) >= 0) & (dist_end_start > dist_end_check));
end
end

function [haveIns, intersection] = check_line_plane_multiPln(start_point, end_point, plane, pt1, pt1o, pt2)
line = end_point - start_point;
t = -(plane(:,4) + dot(repmat(start_point,size(plane,1),1),plane(:,1:3), 2)) ./ dot(repmat(line,size(plane,1),1),plane(:,1:3),2);
intersection = start_point + t .* line;

vector_p1_p2 = pt2 - pt1; % vector AB
vector_p1_pt = intersection - pt1; % vector AM
vector_p1_p1o = pt1o - pt1; % vector AC

AMAB = dot(vector_p1_pt,vector_p1_p2,2);
ABAB = dot(vector_p1_p2,vector_p1_p2,2);
AMAC = dot(vector_p1_pt,vector_p1_p1o,2);
ACAC = dot(vector_p1_p1o,vector_p1_p1o,2);

lie_on_plane = logical((0 < AMAB) & (AMAB < ABAB) & (0 < AMAC) & (AMAC < ACAC));

lies_on_line = check_lie_between_2points(start_point, end_point, intersection);

haveIns = double(logical(lie_on_plane & lies_on_line));
end

