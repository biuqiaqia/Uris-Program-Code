% transform gt_data altitude

clc;
clear;
close all;

load('data\gt_pos\WiFi_sim_2d.mat','gt_data');%static pos file (include correlation positions)

% for WiFi-CP convert 1m above ground to MSL
load('data\model\DTM_HK.mat','DTM_HK');

for idv = 1:size(gt_data,1)
    for idt = 1:size(gt_data{idv,1},1)
         gt_data{idv,1}(idt,4) = Get_HK_MSL(gt_data{idv,1}(idt,2),gt_data{idv,1}(idt,3),DTM_HK)+1;
    end
end



