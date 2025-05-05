clc;
clear;
close all;
D2R = pi/180;
R2D = 180/pi;

%% check data
load('data\sim\static_test_2_sim.mat','data_all','gt_data');

data_all0 = data_all;
vis = [gt_data{1:50,4}];
[~,sort_id]=sortrows(vis');
for idv = 1:1:50
    data_all{1,idv} = data_all0{1,sort_id(idv)};
    data_all{2,idv} = data_all0{2,sort_id(idv)};
    data_all{4,idv} = gt_data{sort_id(idv),4};
    data_all{6,idv} = gt_data{sort_id(idv),3};
end
for idv = 51:1:64
    data_all{4,idv} = gt_data{idv,4};
end

% positioning error estimation
for idv = 1:1:size(data_all,2)
    temp_err = [];
    for idt = 1:1:30
        ls_xyz = llh2xyz([data_all{1,idv}{idt,2}(2),data_all{1,idv}{idt,2}(1),data_all{1,idv}{idt,2}(3)].*[D2R,D2R,1]);
        pos_llh{1,idv}(idt,:) = [data_all{1,idv}{idt,2}(2),data_all{1,idv}{idt,2}(1),data_all{1,idv}{idt,2}(3)];
        gt_xyz = llh2xyz(data_all{2,idv}(idt,1:3).*[D2R,D2R,1]);
        pos_enu = xyz2enu(ls_xyz,gt_xyz);
        temp_err(idt,1) = norm(pos_enu'.*[1,1,0]);
    end
    data_all{5,idv} = rms(temp_err(:,1));
end

f_pos = figure(1);
hold on;
for idv = 1:1:size(data_all,2)
	plot(data_all{2,idv}(1,2),data_all{2,idv}(1,1),'ko','MarkerFaceColor','g','MarkerSize',10);
end
xlabel('Longitude (degree)');
ylabel('Latitude (degree)');
hBase1 = plot_openstreetmap('Alpha',1,'Scale',50000,'BaseUrl', "http://a.tile.openstreetmap.org");
for idv = 1:1:50%size(data_all,2)
    if data_all{4,idv}<=0.15
        plot(data_all{2,idv}(1,2),data_all{2,idv}(1,1),'ko','MarkerFaceColor','g','MarkerSize',10);
%         text(data_all{2,idv}(1,2),data_all{2,idv}(1,1),num2str(idv),'FontSize',10,'Color','k');
%     elseif data_all{4,idv}<=0.3
%         plot(data_all{2,idv}(1,2),data_all{2,idv}(1,1),'ko','MarkerFaceColor','c','MarkerSize',10);
    elseif data_all{4,idv}<=0.53
        plot(data_all{2,idv}(1,2),data_all{2,idv}(1,1),'ko','MarkerFaceColor','b','MarkerSize',10);
    else
        plot(data_all{2,idv}(1,2),data_all{2,idv}(1,1),'ko','MarkerFaceColor','r','MarkerSize',10);
    end
    text(data_all{2,idv}(1,2),data_all{2,idv}(1,1),num2str(idv),'FontSize',15,'Color','k','FontWeight','Bold');
%     text(data_all{2,idv}(1,2),data_all{2,idv}(1,1),num2str(idv),'FontSize',10,'Color','c');
end
% for idv = 51:1:55%size(data_all,2)
%     plot(data_all{2,idv}(1,2),data_all{2,idv}(1,1),'ko','MarkerFaceColor','m','MarkerSize',10);
% %     text(data_all{2,idv}(1,2),data_all{2,idv}(1,1),num2str(idv),'FontSize',15,'Color','k','FontWeight','Bold');
% end
% for idv = 56:1:size(data_all,2)
%     plot(data_all{2,idv}(1,2),data_all{2,idv}(1,1),'ko','MarkerFaceColor','c','MarkerSize',10);
% %     text(data_all{2,idv}(1,2),data_all{2,idv}(1,1),num2str(idv),'FontSize',15,'Color','k','FontWeight','Bold');
% end

figure(2);
hold on;
grid on;
plot([data_all{4,1:50}],'k-','LineWidth',1.2)
ylim([0,1])

%% check specific agent
% ego_id = 1;
% 
% disp(['LS RMSE: ',num2str(data_all{5,ego_id}),', Ang.skymask: ',num2str(data_all{4,ego_id}*100,3)]);
% 
% figure(3);
% hold on;
% plot(pos_llh{ego_id}(:,2),pos_llh{ego_id}(:,1),'r.','MarkerSize',20);
% plot(data_all{2,ego_id}(1,2),data_all{2,ego_id}(1,1),'ko','MarkerSize',10,'MarkerFaceColor','g');
% plot(data_all{2,ego_id}(1,2)+0.0005,data_all{2,ego_id}(1,1)+0.0005,'ko','MarkerSize',10,'MarkerFaceColor','g');
% plot(data_all{2,ego_id}(1,2)-0.0005,data_all{2,ego_id}(1,1)-0.0005,'ko','MarkerSize',10,'MarkerFaceColor','g');
% hBase2 = plot_openstreetmap('Alpha',1,'Scale',50000,'BaseUrl', "http://a.tile.openstreetmap.org");
% plot(pos_llh{ego_id}(:,2),pos_llh{ego_id}(:,1),'r.','MarkerSize',20);
% plot(data_all{2,ego_id}(1,2),data_all{2,ego_id}(1,1),'ko','MarkerSize',10,'MarkerFaceColor','g');
% xlabel('Longitude');
% ylabel('Latitude');
% legend('LS','GT')

% Fig = Skymask_plot_probability_NED([],data_all{6,ego_id},[],2);


