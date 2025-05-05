clc;
clear;
close all;

goGNSS;
para_cl = goGNSS.V_LIGHT/1.023/1000000; 
para_spacing = 1; 

d_phase = [0:2*pi/100:2*pi];

for id = 1:1:100

[dPr_mp(id,1),~] = multipath_pseudorange_delay(0,25,40,30,d_phase(id),para_spacing,para_cl);

end