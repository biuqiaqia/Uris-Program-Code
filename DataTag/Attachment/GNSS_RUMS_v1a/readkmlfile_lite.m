function kml = readkmlfile_lite(kml_path,el_mask,user_xyz)
% function is extracted from Ivan ray-tracing code
D2R = pi/180;
R2D = 180/pi;
S = textread(kml_path,'%s');
idxStart = find(cell2mat(cellfun(@(x) strcmp(x, '<coordinates>'), S, 'UniformOutput', false))) + 1;
idxEnd = find(cell2mat(cellfun(@(x) strcmp(x, '</coordinates>'), S, 'UniformOutput', false))) - 1;
for i = 1:length(idxStart)
    g = cellfun(@(x) strsplit(x, ','), S(idxStart(i):idxEnd(i)), 'UniformOutput',false); 
    g = cell2mat(cellfun(@(x) str2num(x), [g{:}], 'UniformOutput', false));
    g = reshape(g,3,length(g)/3)';
    g = g(:,[2, 1, 3]);
    [x,y,z] = geodetic2ecef(wgs84Ellipsoid,g(:,1),g(:,2),zeros(size(g,1),1));
    e0 = [x,y,z];
    [x,y,z] = geodetic2ecef(wgs84Ellipsoid,g(:,1),g(:,2),g(:,3));
    e = [x,y,z];
    
    p = zeros(size(g,1)-1,4);
    for j = 1:size(g,1)-1
        pt1 = e(j+1,:);
        pt2 = e(j,:);
        pt3 = llh2xyz([g(j,1:2),0.0].*[D2R,D2R,1]);
        [a, b, c, d] = plane_eqn(pt1, pt2, pt3);
        p(j,:) = [a,b,c,d];
    end
    
    kml(i).llh = g;
    kml(i).xyz = e;
    kml(i).xyz0 = e0;
    kml(i).pln = p;
end

% elevation-mask
if ~isempty(el_mask)&&~isempty(user_xyz)
btag = zeros(size(kml,2),1);
for id_b = 1:1:size(kml,2)
	bxyz = mean(kml(id_b).xyz);
	bh = kml(id_b).llh(1,3);
	dist = norm(user_xyz-bxyz);
	bel = asind(bh/dist);
	if bel>el_mask
        btag(id_b,1)=1;
	end
end
kml = kml(btag==1);
end

end

%% plane_eqn
function [a, b, c, d] = plane_eqn(pt1, pt2, pt3)
A = [pt1;pt2;pt3];
A(:,1) = 1;
B = [pt1;pt2;pt3];
B(:,2) = 1;
C = [pt1;pt2;pt3];
C(:,3) = 1;
D = [pt1;pt2;pt3];
D = D * -1;

a = det(A);
b = det(B);
c = det(C);
d = det(D);
end