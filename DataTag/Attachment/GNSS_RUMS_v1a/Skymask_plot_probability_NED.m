function Fig = Skymask_plot_probability_NED(SV_elaz,skymask,outpath,epoch)
% -------------------------------------------------------
% Function for plotting skymask with satellite visibility
% In NED coordinate frame
% Figure also saved in the direction
%
% by GH.Zhang 2018/12/14
% guo-hao.zhang@connect.polyu.hk
% -------------------------------------------------------
% Input:    SV_elaz(column 1)	-   Satellite PRN number
%           SV_elaz(column 2)	-   Satellite elevation
%           SV_elaz(column 3)	-   Satellite azimuth in NED frame
%           SV_elaz(column 4)	-   Satellite visibility (0-NLOS/1-LOS)
%           skymask             -   Building boundary elevation in NED frame
%           outpath             -   Output file direction
%           epoch               -   Current epoch (just to use as figure title)
%
% Output:   Fig     -   Skyplot figure

% construct colormap
cmap(:,1) = (0:1:240)'/360;
cmap(:,2:3) = ones(241,2);
for idco =1:1:241
    cmap2(idco,:)=hsv2rgb(cmap(242-idco,:));
end 

ref_data = [];
for az = 1:1:360
    mask_el = skymask(az,1);
    for el = 0:1:90
        if el <= mask_el
            tag = 2;
        else
            tag = 1;
        end
        ref_data = [ref_data;el,az,tag];
    end
end
elevation_ref = ref_data(:,1);
azimuth_ref = ref_data(:,2);
signaltype_ref = ref_data(:,3);

Fig = figure(epoch);
a = 1;
patch([-a -a a  a],...
      [-a  a a -a],[239/255 239/255 239/255],'EdgeColor','none');

radius = 0.9;
LOSidx = find(signaltype_ref==1);
NLOSidx = find(signaltype_ref==2 | signaltype_ref==3);
hold on;
h=plot(radius*(ones(size(LOSidx,1),1)-elevation_ref(LOSidx)./90).*sind(azimuth_ref(LOSidx)),radius*(ones(size(LOSidx,1),1)-elevation_ref(LOSidx)./90).*cosd(azimuth_ref(LOSidx)),'.','MarkerEdgeColor',[0.99 0.99 0.99]);
    plot(radius*(ones(size(NLOSidx,1),1)-elevation_ref(NLOSidx)./90).*sind(azimuth_ref(NLOSidx)),radius*(ones(size(NLOSidx,1),1)-elevation_ref(NLOSidx)./90).*cosd(azimuth_ref(NLOSidx)),'.','MarkerEdgeColor',[0.5 0.5 0.5]);

idflg=1;
figflg=1;
flag=1;

hlinex=[-0.9 0.9]; hliney=[0 0]; vlinex=[0 0]; vliney=[-0.9 0.9];
arg=[0:100]*2*pi/100;
% if figflg
%close
% end
plot(0.9*sin(arg),0.9*cos(arg),'k--',...
     0.8*sin(arg),0.8*cos(arg),'r-',...
     0.7*sin(arg),0.7*cos(arg),'k--',...
     0.6*sin(arg),0.6*cos(arg),'k--',...
     0.5*sin(arg),0.5*cos(arg),'k--',...
     0.4*sin(arg),0.4*cos(arg),'k--',...
     0.3*sin(arg),0.3*cos(arg),'k--',...
     0.2*sin(arg),0.2*cos(arg),'k--',...
     0.1*sin(arg),0.1*cos(arg),'k--',...
     hlinex,hliney,'k-',vlinex,vliney,'k-')
% plot(0.9*sin(arg),0.9*cos(arg),'k--',...
%      0.6*sin(arg),0.6*cos(arg),'k--',...
%      0.3*sin(arg),0.3*cos(arg),'k--',...
%      hlinex,hliney,'k-',vlinex,vliney,'k-')
axis('square')
axis('off')
% text(-0.1,0.85,'N','FontSize',14)
% if ~isempty(epoch)
%     title(['Skymask - ',num2str(epoch)]);
% end

% --- plot nSV position ---
if isempty(SV_elaz)
    Fig = nan;
elseif size(SV_elaz,2)==1
    for prn = 1:1:size(SV_elaz,1)
        if ~isempty(SV_elaz{prn,1})
            if prn<=32
                prn_text = ['G',num2str(prn)];
            elseif prn>=87
                prn_text = ['B',num2str(prn-86)];
            end
            plot(radius*(-SV_elaz{prn,1}(:,1)./90+1).*sind(SV_elaz{prn,1}(:,2)),...
                 radius*(-SV_elaz{prn,1}(:,1)./90+1).*cosd(SV_elaz{prn,1}(:,2)),...
                 '.','MarkerEdgeColor',[0,0,1],'MarkerSize',20);  
            text(radius*(1-SV_elaz{prn,1}(1,1)/90).*sind(SV_elaz{prn,1}(1,2))+0.05,...
                 radius*(1-SV_elaz{prn,1}(1,1)/90).*cosd(SV_elaz{prn,1}(1,2)),...
                 prn_text,'FontSize',15,'Color',[0,0,0],'FontWeight','Bold');
        end
    end
else
%     max_score = 1;
%     for idsvp = 1:1:size(SV_elaz,1)
%         if SV_elaz(idsvp,4)>1.01
%             return
%         elseif SV_elaz(idsvp,4)>1
%             SV_elaz(idsvp,4)=1;
%         end
%         if isnan(SV_elaz(idsvp,4))
%             plot(radius*(1-SV_elaz(idsvp,2)/90).*sind(SV_elaz(idsvp,3)),...
%                  radius*(1-SV_elaz(idsvp,2)/90).*cosd(SV_elaz(idsvp,3)),...
%                  'ko','MarkerSize',20,'LineWidth',2); 
%         else
%             plot(radius*(1-SV_elaz(idsvp,2)/90).*sind(SV_elaz(idsvp,3)),...
%                  radius*(1-SV_elaz(idsvp,2)/90).*cosd(SV_elaz(idsvp,3)),...
%                  '.','MarkerEdgeColor',...
%                  hsv2rgb((240-interp1(0:max_score/(length(0:240)-1):max_score,0:240,1-SV_elaz(idsvp,4)))/360,1,1),...
%                  'MarkerSize',75);  
%         end
%     end
    for idsvp = 1:1:size(SV_elaz,1)
        if SV_elaz(idsvp,4) == 1
            plot(radius*(1-SV_elaz(idsvp,2)/90).*sind(SV_elaz(idsvp,3)),...
                 radius*(1-SV_elaz(idsvp,2)/90).*cosd(SV_elaz(idsvp,3)),...
                 'ko','MarkerFaceColor',[0,1,0],'MarkerSize',20);  
        elseif SV_elaz(idsvp,4) == 0
            plot(radius*(1-SV_elaz(idsvp,2)/90).*sind(SV_elaz(idsvp,3)),...
                 radius*(1-SV_elaz(idsvp,2)/90).*cosd(SV_elaz(idsvp,3)),...
                 'ko','MarkerFaceColor',[1,0,0],'MarkerSize',20);
        end
    end
    for idsvp = 1:1:size(SV_elaz,1)
        if SV_elaz(idsvp,1)>=10
            text(radius*(1-SV_elaz(idsvp,2)/90).*sind(SV_elaz(idsvp,3))-0.08,...
                 radius*(1-SV_elaz(idsvp,2)/90).*cosd(SV_elaz(idsvp,3)),...
                 num2str(SV_elaz(idsvp,1)),'FontSize',15,'Color',[0,0,0],'FontWeight','Bold');
        else
            text(radius*(1-SV_elaz(idsvp,2)/90).*sind(SV_elaz(idsvp,3))-0.04,...
                 radius*(1-SV_elaz(idsvp,2)/90).*cosd(SV_elaz(idsvp,3)),...
                 num2str(SV_elaz(idsvp,1)),'FontSize',15,'Color',[0,0,0],'FontWeight','Bold');
        end
    end
%     colormap(cmap2);
%     caxis([0 max_score]);
%     colorbar;
end

if ~isempty(outpath)
    saveas(Fig,[outpath,num2str(epoch)],'jpeg');
%     close(Fig)
end


