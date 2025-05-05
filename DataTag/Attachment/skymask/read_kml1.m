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
        disp(i);
        for j = 1:size(C,1)
            xyz(j,:) = llh2xyz([C(j,1:2)*D2R C(j,3)]);
        end
        building(i,1).v_llh = C;
        building(i,1).v_xyz = xyz;
        building(i,1).center = mean(C);
        building(i,1).center_xyz = mean(xyz);
        clear xyz, C;
    end
end

