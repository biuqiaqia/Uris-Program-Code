1. Generate skymask
1.1 Run 'main_gen_skymask.m'
1.2 Un-comment kml block you want
    Block example
	kml_file = 'kml\Barcelona_P9_P10.kml'; % kml file
	pos = [41.396205	2.153784]; % ground truth position
	d = 220; % 
	lat2m = distance(pos(1),pos(2),pos(1)+1,pos(2),wgs84Ellipsoid); m2lat = 1/lat2m;
	lon2m = distance(pos(1),pos(2),pos(1),pos(2)+1,wgs84Ellipsoid); m2lon = 1/lon2m;
	max_lat = pos(1)+d*m2lat; % calculate upper boundary
	min_lat = pos(1)-d*m2lat; % calculate lower boundary
	max_lon = pos(2)+d*m2lon; % calculate right boundary
	min_lon = pos(2)-d*m2lon; % calculate leftboundary
1.3 If temp file found, it will load itself
1.4 Wait until finish generating.


2. Merge formatted skymask
2.1 Run 'skymask_merge.m'
2.2 Un-comment kml line
    Kml line example:
	building_model_path = 'kml\Frankfurt_P12.kml';
2.3 It will find all temp file and merge to final skymask saved in 'mat\' directory.