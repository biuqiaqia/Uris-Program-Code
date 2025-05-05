function CN0_open = opensky_CN0_modeling(PRN,fit_para,el,type)

% calculating opensky C/N0 by data modeling
% type 1-linear / 2-logarithm

if PRN<=32
    cnr = polyval(fit_para.gps,el);
elseif ismember(PRN,fit_para.bds_type.geo)
    cnr = polyval(fit_para.bds_geo,el);
elseif ismember(PRN,fit_para.bds_type.igso)
    cnr = polyval(fit_para.bds_igso,el);
elseif ismember(PRN,fit_para.bds_type.meo)
    cnr = polyval(fit_para.bds_meo,el);
end

if type == 1
    CN0_open = cnr;
elseif type == 2
    CN0_open = 10*log10(cnr);
end

