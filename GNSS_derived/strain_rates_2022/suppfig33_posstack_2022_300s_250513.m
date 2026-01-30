%% Stevens et al. [2015] Supplementary Info FIGURE 6 --> 2022 displacements
% GPS vs NIF: all stations data in 3 long columns
% 11 June 2021: change for flowline stations in 2012 drainage, remove NIF
% outputs
% 13 Oct 2021: remove 3\sigma outliers for the rough day boundary
% 2023 July 17: LAS revisions for reviews
% 2023 Oct 10: LAS change for 2022 lake6 array; TRACK positions

clear all; close all
%% load 2022 meltseason positions
load A2022_unfixedBF2_atm_datarejected_260112.mat
station_names = {A2022.name}; sites = station_names;
num_sta = length(station_names)-2; % number of stations; omit SQ33
station_names_flip = vertcat(station_names(22),station_names(21),station_names(20),station_names(19),...
    station_names(18),station_names(17),station_names(16),...
    station_names(15),station_names(14),station_names(13),...
    station_names(12),station_names(11),station_names(10),...
    station_names(9),station_names(8),station_names(7),...
    station_names(6),station_names(5),station_names(4),...
    station_names(3),station_names(2),station_names(1));
save station_names.mat station_names
station_names = upper(station_names);
open_time_L1A_2022 = 195.0972;
load B2022R_BF2_UP4_sZEROm_OUT6hr_clean_260113.mat % cleaned up positions

%% positions into along- and across-flow directions based on springtime vels
springtime = [150, 160]; 
interptime_15sec = 1.736111100001381e-04; % 15-sec in decimal day
ti22_15sec = (springtime(1):interptime_15sec:springtime(2))'; % 15-second time vector

% flow angle of springtime velocity 
for i=1:1:num_sta 
    % time indices
    idx_time(i,1) = find((springtime(1)-B2022(i).neu_iu(:,1))<=0,1); % start
    idx_time(i,2) = find((springtime(2)-B2022(i).neu_iu(:,1))<=0,1); % stop

    % calculate velocity over DOY 150-160
    % horizontal vel: quick distance moved method
    bg_moved(i,1) = sqrt(((nanmean(B2022(i).neu_iu(idx_time(i,2),2))-nanmean(B2022(i).neu_iu(idx_time(i,1),2))).^2) +...
        ((nanmean(B2022(i).neu_iu(idx_time(i,2),3))-nanmean(B2022(i).neu_iu(idx_time(i,1),3))).^2)); % [ m ]
    % velocity 
    bg_moved(i,2) = bg_moved(i,1)/(springtime(2)-springtime(1)); % [ m/day ]
    bg_moved(i,3) = (bg_moved(i,1)/(springtime(2)-springtime(1)))*365.25; % [ m/yr ]
    % flow direction
    bg_moved(i,4) = (atan2d((B2022(i).neu_iu(idx_time(i,2),2)-B2022(i).neu_iu(idx_time(i,1),2)),...
        (B2022(i).neu_iu(idx_time(i,2),3)-B2022(i).neu_iu(idx_time(i,1),3)))) ;  % 0 degrees = east !!!

    % horizontal vel: linear fit to DOY 150-160 method
    B2022(i).bg_displ(:,1) = sqrt(((B2022(i).neu_iu(idx_time(i,1):idx_time(i,2),2)-B2022(i).neu_iu(idx_time(i,1),2)).^2) + ...
        ((B2022(i).neu_iu(idx_time(i,1):idx_time(i,2),3)-B2022(i).neu_iu(idx_time(i,1),3)).^2)); % [ m ]
    lin_fit(i,:) = polyfit(B2022(i).neu_iu(idx_time(i,1):idx_time(i,2),1)-ti22_15sec(1),...
        B2022(i).bg_displ(:,1), 1); 
    % velocity/flow direction -- linear fit method
    bg_moved(i,5) = B2022(i).bg_displ(end,1); % [ m ]
    bg_moved(i,6) = lin_fit(i,1); % [ m/day ]
    bg_moved(i,7) = lin_fit(i,1)*365.25; % [ m/yr ]
    if bg_moved(i,4) >= 0
        bg_moved(i,8) = 360 - bg_moved(i,4) + 90 ; % 0 degrees = north !!!
    else if bg_moved(i,4) < 0 
        bg_moved(i,8) = (-1.*bg_moved(i,4)) + 90 ; % 0 degrees = north !!!
    end
    end

    % vertical vel: linear fit to DOY 150-160 method
    B2022(i).bg_displ(:,2) = B2022(i).neu_iu(idx_time(i,1):idx_time(i,2),4)-B2022(i).neu_iu(idx_time(i,1),4); % [ m ]
    lin_fit_vertical(i,:) = polyfit(B2022(i).neu_iu(idx_time(i,1):idx_time(i,2),1)-ti22_15sec(1), B2022(i).bg_displ(:,2), 1); 
    % velocity/flow direction -- linear fit method
    bg_moved(i,9) = B2022(i).bg_displ(end,2); % [ m ]
    bg_moved(i,10) = lin_fit_vertical(i,1); % [ m/day ]
    bg_moved(i,11) = lin_fit_vertical(i,1)*365.25; % [ m/yr ]
end

%% horizontal positions into along- and across-flow directions
for i=1:1:num_sta 
   % flow angle 
   flow_angle(i,1) =  bg_moved(i,8); % flow direction [ deg ]  [0=north]
   flow_angle(i,2) =  (360-bg_moved(i,8)) + 90; % flow direction [ deg ]  [0=east]
   % displacement from t_0
   B2022(i).neu_displ(:,1) = B2022(i).neu_iu(:,1);
   B2022(i).neu_displ(:,2) = B2022(i).neu_iu(:,2) - B2022(i).neu_iu(1,2); % N
   B2022(i).neu_displ(:,3) = B2022(i).neu_iu(:,3) - B2022(i).neu_iu(1,3); % E
   B2022(i).neu_displ(:,4) = B2022(i).neu_iu(:,4) - B2022(i).neu_iu(1,4); % U
   % displacement from t_0 w/flowangle
   B2022(i).neu_displ_flow(:,1) = B2022(i).neu_iu(:,1);
   B2022(i).neu_displ_flow(:,2) = -1.*((B2022(i).neu_iu(:,2) - B2022(i).neu_iu(1,2)).*cosd(flow_angle(i,2))) + ...
       ((B2022(i).neu_iu(:,3) - B2022(i).neu_iu(1,3)).*sind(flow_angle(i,2))); % ACROSS
   B2022(i).neu_displ_flow(:,3) = ((B2022(i).neu_iu(:,3) - B2022(i).neu_iu(1,3)).*cosd(flow_angle(i,2))) + ...
       ((B2022(i).neu_iu(:,2) - B2022(i).neu_iu(1,2)).*sind(flow_angle(i,2))); % ALONG
   B2022(i).neu_displ_flow(:,4) = B2022(i).neu_iu(:,4) - B2022(i).neu_iu(1,4); % UP  
end
% check fig
% figure; clf;
% plot(B2022(14).neu_displ_flow(:,1),B2022(14).neu_displ_flow(:,2),'.'); 
% 
% SQ35 gnarly edge
for i=19
time_range = [209.01 209.2]; % [decimal days]
    [II] = find(B2022(i).neu_displ_flow(2:end,1)>(time_range(1)) & ...
           B2022(i).neu_displ_flow(2:end,1)<(time_range(2)));    
    B2022(i).neu_displ_flow(II,1:4) = NaN;
end
FLAG = isnan(B2022(19).neu_displ_flow(:,1));
B2022(19).neu_displ_flow(FLAG,:) = [];

%% load north lake geographic files and morlighem bed relative to North Lake
origin = [68.72, -49.53]; % M1 moulin
load polarstereo_stations_2022_short.mat
% % polarstereo conversion needed values
radius=6378137.0;    eccen=0.08181919;    lat_true=70;    lon_posy=-45; 
[moulin_x, moulin_y]=polarstereo_fwd(origin(1,1), origin(1,2),radius,eccen,lat_true,lon_posy);
% load bedmap
load('BMv5_for_nevis_catchment.mat'); % origin = M1 moulin 
load cmapland2.mat % bed topo colormap
% load lake boundaries, centre points, drainage dates from FASTER for 2022
load('environs_lakes_2022B_250416.mat') % drainage dates & boundaries
environs_lakes_2022 = environs_lakes; environs_lakes_2022B = environs_lakes; 

% dock at eel pond :)
sky_blue = [111, 169, 228]./255; % SNL
metal = [87, 115, 131]./255; % NNL
oar = [251, 219, 154]./255; % NENL
handle = [161, 37, 49]./255; % NL
dark_oar = [164, 114, 63]./255;
% more colors
ice_sheet_contours = [0.7 0.7 0.7];
plum = [211, 160, 211, 100]./255;
lake_blue = [0    0.4470    0.7410   0.08]; % fourth value is FaceAlpha built in!
goldenrod = [1.0 0.8 0.0];
icy_plum = [203, 183, 227]./255;

%% scrap figure for colors
Fig1 = figure(1); clf; 
set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',[5 2 19 14]);
axe1 = axes('Position',[0.075 0.70 0.9 0.25],'Box','on','NextPlot','add','XTickLabels',[]);
axes(axe1)
h1=plot([168.85 168.85], [0 100],'-'); c1 = get(h1,'Color'); hold on;
h2=plot([169.85 168.85], [0 100],'-'); c2 = get(h2,'Color');
h3=plot([170.85 168.85], [0 100],'-'); c3 = get(h3,'Color');
h4=plot([171.85 168.85], [0 100],'-'); c4 = get(h4,'Color');
h5=plot([172.85 168.85], [0 100],'-'); c5 = get(h5,'Color');
h6=plot([173.85 168.85], [0 100],'-'); c6 = get(h6,'Color');
h7=plot([174.85 168.85], [0 100],'-'); c7 = get(h7,'Color');
h8=plot([175.85 168.85], [0 100],'-'); c8 = get(h8,'Color');
h9=plot([176.85 168.85], [0 100],'-'); c9 = get(h9,'Color');

%% Figure 300s 2022
figure(2); clf;
set(gcf,'Units','centimeters','Position',1.5.*[0.5 0 19.3 20]);
    
xmin = 209.3; xmax = 211.0; 
ymin = -0.10; ymax = 1.75;
offset_all = 5.25;
    
    axe0 = axes('Position',[0.07+0.29+0.045 0.7325 0.225 0.2125],'Box','on','NextPlot','add');

    axe1 = axes('Position',[0.059 0.04 0.29 0.94],'Box','on','NextPlot','add');
    ylabel('GNSS Station','FontSize',11);
    set(gca, 'XTickLabelMode', 'auto');
    set(gca,'YTick',0:0.25:ymax,'YTickLabel',{'SQ37';'SQ36';'SQ35';'SQ34';'SQ31';...
        '';'';'';'';'';'';'';'';'';''})
    
    axe2 = axes('Position',[0.055+0.30 0.04 0.29 0.94],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'auto');
    set(gca,'YTick',0:0.25:ymax,'YTickLabel',[]);
    xlabel('Day of Year, 2022  [ UTC ]','FontSize',11,'Fontname','Avenir')

    axe3 = axes('Position',[0.053+0.30+0.30 0.04 0.29 0.94],'Box','on','NextPlot','add');
    ylabel('Displacement  [ m ]','FontSize',11,'Fontname','Avenir');
    set(gca, 'XTickLabelMode', 'auto');
    set(gca,'yaxislocation','right','ytick',0:0.25:ymax)
    
    TriangleSize = 7;
    offset = 0.25; % displacements vertical offset
    offset_v = offset.*(0:1:21); % displacements y-intercept

% begin GPS displacements
for i=17:18 % 1350s num_sta
    % find within the hour of xmin (over a 1-hr window)
    time_before = 6/24; % [decimal days]
    [II] = find(B2022(i).neu_displ_flow(2:end,1)>(xmin-(time_before/2)) & ...
          B2022(i).neu_displ_flow(2:end,1)<(xmin+(time_before/2)));    
    pos_30min_before_open_time_L1A_2022(i,:) = nanmean(B2022(i).neu_displ_flow(II,1:4));
    % find background velocity 1 day before xmin (over a 24-hr window)
    time_before = 1; % [decimal days]
    [I] = find(B2022(i).neu_displ_flow(2:end,1)>(xmin-(time_before)) & ...
          B2022(i).neu_displ_flow(2:end,1)<(xmin+0.2));
    vel_30min_before(i,1) = B2022(i).neu_displ_flow(I(1),1); % time
        % slopes
        s_across = polyfit(B2022(i).neu_displ_flow(I,1),B2022(i).neu_displ_flow(I,2),1); % slope in across [m/day]
        s_along = polyfit(B2022(i).neu_displ_flow(I,1),B2022(i).neu_displ_flow(I,3),1); % slope in along [m/day]
        s_up = polyfit(B2022(i).neu_displ_flow(I,1),B2022(i).neu_displ_flow(I,4),1); % slope in up [m/day]
        % velocities
        vel_30min_before(i,2) = s_along(1); % [ m/d ]
        vel_30min_before(i,3) = s_across(1); % [ m/d ]
        vel_30min_before(i,4) = s_up(1); % [ m/d ]

  axes(axe1);% E
  title('Flowline (276-280˚) displacement [ m ]','Fontname','Avenir','FontSize',11);
  plot_data_along = offset_all + (B2022(i).neu_displ_flow(:,3)-pos_30min_before_open_time_L1A_2022(i,3)) - ...
                    offset_v(i+1) - (vel_30min_before(i,2).*(B2022(i).neu_displ_flow(:,1)-xmin));
  plot(B2022(i).neu_displ_flow(:,1),plot_data_along,'.','markerSize',2.5); 
  hold on; grid on; set(gca,'XMinorGrid','on','YMinorGrid','on','LineWidth',1.1,'FontSize',10,'FontName','Avenir')
  %text(xmin+0.03, ymax-0.07, 'a','FontWeight','bold','FontSize',12);
  xlim([xmin xmax]); ylim([ymin ymax])
  box on

  axes(axe2) %N
  title('Across-flow (006-010˚) displacement [ m ]','Fontname','Avenir','FontSize',11)
  plot_data_across = offset_all + (B2022(i).neu_displ_flow(:,2)-pos_30min_before_open_time_L1A_2022(i,2)) - ...
      offset_v(i+1) - (vel_30min_before(i,3).*(B2022(i).neu_displ_flow(:,1)-xmin));
  plot(B2022(i).neu_displ_flow(:,1),plot_data_across,'.','markerSize',2.5); 
  hold on; grid on; set(gca,'XMinorGrid','on','YMinorGrid','on','LineWidth',1.1,'FontSize',10,'FontName','Avenir')
  %text(xmin+0.03, ymax-0.07, 'b','FontWeight','bold','FontSize',12);
  xlim([xmin xmax]); ylim([ymin ymax])
  box on
  
  axes(axe3) %U
  title('Vertical displacement [ m ]','Fontname','Avenir','FontSize',11)
  plot_data_up = offset_all + (B2022(i).neu_displ_flow(:,4)-pos_30min_before_open_time_L1A_2022(i,4)) - ...
      offset_v(i+1) - (vel_30min_before(i,4).*(B2022(i).neu_displ_flow(:,1)-xmin));
  plot(B2022(i).neu_displ_flow(:,1),plot_data_up,'.','markerSize',2.5); 
  hold on; grid on; set(gca,'XMinorGrid','on','YMinorGrid','on','LineWidth',1.1,'FontSize',10,'FontName','Avenir')
  %text(xmin+0.03, ymax-0.07, 'c','FontWeight','bold','FontSize',12);
  xlim([xmin xmax]); ylim([ymin ymax])
  set(gca,'ytick',0:0.25:1.25)
end

% begin GPS displacements
for i=19 % SQ35 (avoiding data gap at xmin)
    % find within the hour of xmin (over a 1-hr window)
    time_before = 8/24; % [decimal days]
    [II] = find(B2022(i).neu_displ_flow(2:end,1)>(208.5) & ...
          B2022(i).neu_displ_flow(2:end,1)<(208.9));    
    pos_30min_before_open_time_L1A_2022(i,:) = nanmean(B2022(i).neu_displ_flow(II,1:4));
    % find background velocity 1 day before xmin (over a 24-hr window)
    time_before = 1; % [decimal days]
    [I] = find(B2022(i).neu_displ_flow(2:end,1)>(208.0) & ...
          B2022(i).neu_displ_flow(2:end,1)<(208.9));
    vel_30min_before(i,1) = B2022(i).neu_displ_flow(I(1),1); % time
        % slopes
        s_across = polyfit(B2022(i).neu_displ_flow(I,1),B2022(i).neu_displ_flow(I,2),1); % slope in across [m/day]
        s_along = polyfit(B2022(i).neu_displ_flow(I,1),B2022(i).neu_displ_flow(I,3),1); % slope in along [m/day]
        s_up = polyfit(B2022(i).neu_displ_flow(I,1),B2022(i).neu_displ_flow(I,4),1); % slope in up [m/day]
        % velocities
        vel_30min_before(i,2) = s_along(1); % [ m/d ]
        vel_30min_before(i,3) = s_across(1); % [ m/d ]
        vel_30min_before(i,4) = s_up(1); % [ m/d ]

  axes(axe1);% E
  title('Flowline (276-280˚) displacement [ m ]','Fontname','Avenir','FontSize',11);
  plot_data_along = offset_all + (B2022(i).neu_displ_flow(:,3)-pos_30min_before_open_time_L1A_2022(i,3)) - ...
                    offset_v(i+2) - (vel_30min_before(i,2).*(B2022(i).neu_displ_flow(:,1)-xmin));
  plot(B2022(i).neu_displ_flow(:,1),plot_data_along,'.','markerSize',2.5); 
  hold on; grid on; set(gca,'XMinorGrid','on','YMinorGrid','on','LineWidth',1.1,'FontSize',10,'FontName','Avenir')
  xlim([xmin xmax]); ylim([ymin ymax])
  box on

  axes(axe2) %N
  title('Across-flow (006-010˚) displacement [ m ]','Fontname','Avenir','FontSize',11)
  plot_data_across = offset_all + (B2022(i).neu_displ_flow(:,2)-pos_30min_before_open_time_L1A_2022(i,2)) - ...
      offset_v(i+1) - (vel_30min_before(i,3).*(B2022(i).neu_displ_flow(:,1)-xmin));
  plot(B2022(i).neu_displ_flow(:,1),plot_data_across,'.','markerSize',2.5); 
  hold on; grid on; set(gca,'XMinorGrid','on','YMinorGrid','on','LineWidth',1.1,'FontSize',10,'FontName','Avenir')
  xlim([xmin xmax]); ylim([ymin ymax])
  box on
  
  axes(axe3) %U
  title('Vertical displacement [ m ]','Fontname','Avenir','FontSize',11)
  plot_data_up = offset_all + (B2022(i).neu_displ_flow(:,4)-pos_30min_before_open_time_L1A_2022(i,4)) - ...
      offset_v(i+1) - (vel_30min_before(i,4).*(B2022(i).neu_displ_flow(:,1)-xmin));
  plot(B2022(i).neu_displ_flow(:,1),plot_data_up,'.','markerSize',2.5); 
  hold on; grid on; set(gca,'XMinorGrid','on','YMinorGrid','on','LineWidth',1.1,'FontSize',10,'FontName','Avenir')
  xlim([xmin xmax]); ylim([ymin ymax])
  set(gca,'ytick',0:0.25:1.25)
end

% begin GPS displacements
for i=20:21 % 1350s num_sta
    % find within the hour of xmin (over a 1-hr window)
    time_before = 1/24; % [decimal days]
    [II] = find(B2022(i).neu_displ_flow(2:end,1)>(xmin-(time_before/2)) & ...
          B2022(i).neu_displ_flow(2:end,1)<(xmin+(time_before/2)));    
    pos_30min_before_open_time_L1A_2022(i,:) = nanmean(B2022(i).neu_displ_flow(II,1:4));
    % find background velocity 1 day before xmin (over a 24-hr window)
    time_before = 1.0; % [decimal days]
    [I] = find(B2022(i).neu_displ_flow(2:end,1)>(xmin-(time_before)) & ...
          B2022(i).neu_displ_flow(2:end,1)<(xmin+0.2));
    vel_30min_before(i,1) = B2022(i).neu_displ_flow(I(1),1); % time
        % slopes
        s_across = polyfit(B2022(i).neu_displ_flow(I,1),B2022(i).neu_displ_flow(I,2),1); % slope in across [m/day]
        s_along = polyfit(B2022(i).neu_displ_flow(I,1),B2022(i).neu_displ_flow(I,3),1); % slope in along [m/day]
        s_up = polyfit(B2022(i).neu_displ_flow(I,1),B2022(i).neu_displ_flow(I,4),1); % slope in up [m/day]
        % velocities
        vel_30min_before(i,2) = s_along(1); % [ m/d ]
        vel_30min_before(i,3) = s_across(1); % [ m/d ]
        vel_30min_before(i,4) = s_up(1); % [ m/d ]

  axes(axe1);% E
  title('Flowline (276-280˚) displacement [ m ]','Fontname','Avenir','FontSize',11);
  plot_data_along = offset_all + (B2022(i).neu_displ_flow(:,3)-pos_30min_before_open_time_L1A_2022(i,3)) - ...
                    offset_v(i+1) - (vel_30min_before(i,2).*(B2022(i).neu_displ_flow(:,1)-xmin));
  plot(B2022(i).neu_displ_flow(:,1),plot_data_along,'.','markerSize',2.5); 
  hold on; grid on; set(gca,'XMinorGrid','on','YMinorGrid','on','LineWidth',1.1,'FontSize',10,'FontName','Avenir')
  text(xmin+0.03, ymax-0.03, 'a','FontWeight','bold','FontSize',12);
  xlim([xmin xmax]); ylim([ymin ymax])
  box on

  axes(axe2) %N
  title('Across-flow (006-010˚) displacement [ m ]','Fontname','Avenir','FontSize',11)
  plot_data_across = offset_all + (B2022(i).neu_displ_flow(:,2)-pos_30min_before_open_time_L1A_2022(i,2)) - ...
      offset_v(i+1) - (vel_30min_before(i,3).*(B2022(i).neu_displ_flow(:,1)-xmin));
  plot(B2022(i).neu_displ_flow(:,1),plot_data_across,'.','markerSize',2.5); 
  hold on; grid on; set(gca,'XMinorGrid','on','YMinorGrid','on','LineWidth',1.1,'FontSize',10,'FontName','Avenir')
  text(xmin+0.03, ymax-0.03, 'b','FontWeight','bold','FontSize',12);
  xlim([xmin xmax]); ylim([ymin ymax])
  box on

  axes(axe3) %U
  title('Vertical displacement [ m ]','Fontname','Avenir','FontSize',11)
  plot_data_up = offset_all + (B2022(i).neu_displ_flow(:,4)-pos_30min_before_open_time_L1A_2022(i,4)) - ...
      offset_v(i+1) - (vel_30min_before(i,4).*(B2022(i).neu_displ_flow(:,1)-xmin));
  plot(B2022(i).neu_displ_flow(:,1),plot_data_up,'.','markerSize',2.5); 
  hold on; grid on; set(gca,'XMinorGrid','on','YMinorGrid','on','LineWidth',1.1,'FontSize',10,'FontName','Avenir')
  text(xmin+0.03, ymax-0.03, 'c','FontWeight','bold','FontSize',12);
  xlim([xmin xmax]); ylim([ymin ymax])
  set(gca,'ytick',0:0.25:1.25)
end
    
axes(axe0);
hold on;
text(36.25, -21.5, 'd','FontWeight','bold','FontSize',11);
% lake names
text(39, -24.25, 'L3A','FontSize',11,'FontName','AvenirBold','FontAngle','italic','Color','w','FontWeight','bold')
text(37.25, -28.0, 'L3B','FontSize',11,'FontName','Avenir','FontAngle','italic')
% FASTER lakes
% lakes: FASTER S2 2022
[HF_lakes_to_plot_22] = find(environs_lakes_2022.laketypeing_dates(:,7)==1); % just HF lakes   
[MU_lakes_to_plot_22] = find(environs_lakes_2022.laketypeing_dates(:,7)==2); % just moulin-draining lakes 
[lakes_to_plot_22] = find(~isnan(environs_lakes_2022.laketypeing_dates(:,7))); % all good lakes   
[NE_lakes_to_plot_22] = find(environs_lakes_2022.laketypeing_dates(:,7)==5); % just no-exit freezers
% all lakes in blue
for i=1:1:length(lakes_to_plot_22)
    fill3(environs_lakes_2022B.boundaries(lakes_to_plot_22(i)).XY_km_local(:,1), ...
        environs_lakes_2022B.boundaries(lakes_to_plot_22(i)).XY_km_local(:,2), ...
        -3000.*ones(length(environs_lakes_2022B.boundaries(lakes_to_plot_22(i)).XY_km_local(:,2)),1),... 
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',lake_blue(1:3));
end
% HF events
for i=1:1:length(HF_lakes_to_plot_22)
    fill3(environs_lakes_2022B.boundaries(HF_lakes_to_plot_22(i)).XY_km_local(:,1), ...
        environs_lakes_2022B.boundaries(HF_lakes_to_plot_22(i)).XY_km_local(:,2), ...
        -3000.*ones(length(environs_lakes_2022B.boundaries(HF_lakes_to_plot_22(i)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',goldenrod);
end
% MU events
for i=1:1:length(MU_lakes_to_plot_22)
    fill3(environs_lakes_2022B.boundaries(MU_lakes_to_plot_22(i)).XY_km_local(:,1), ...
        environs_lakes_2022B.boundaries(MU_lakes_to_plot_22(i)).XY_km_local(:,2), ...
         -3000.*ones(length(environs_lakes_2022B.boundaries(MU_lakes_to_plot_22(i)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',[0.4 0.8 0.8]);
end
% No-exit freezers 
for i=1:1:length(NE_lakes_to_plot_22)
    fill3(environs_lakes_2022B.boundaries(NE_lakes_to_plot_22(i)).XY_km_local(:,1), ...
        environs_lakes_2022B.boundaries(NE_lakes_to_plot_22(i)).XY_km_local(:,2), ...
        -3000.*ones(length(environs_lakes_2022B.boundaries(NE_lakes_to_plot_22(i)).XY_km_local(:,2)),1),...
        'o','MarkerSize',1,'EdgeColor','none','FaceColor',plum(1:3));
end
% GNSS stations
text(xy_sta_22_short(17,1)-0.40,xy_sta_22_short(17,2)-0.5,station_names(17),'FontSize',9,'Rotation',0);
text(xy_sta_22_short(18:21,1)-0.40,xy_sta_22_short(18:21,2)-0.5,station_names(20:23),'FontSize',9,'Rotation',0);
plot(xy_sta_22_short(17,1),xy_sta_22_short(17,2),'k^','MarkerSize',TriangleSize-0,'MarkerFaceColor',c1);
plot(xy_sta_22_short(18,1),xy_sta_22_short(18,2),'k^','MarkerSize',TriangleSize-0,'MarkerFaceColor',c2);
plot(xy_sta_22_short(19,1),xy_sta_22_short(19,2),'k^','MarkerSize',TriangleSize-0,'MarkerFaceColor',c3);
plot(xy_sta_22_short(20,1),xy_sta_22_short(20,2),'k^','MarkerSize',TriangleSize-0,'MarkerFaceColor',c4);
plot(xy_sta_22_short(21,1),xy_sta_22_short(21,2),'k^','MarkerSize',TriangleSize-0,'MarkerFaceColor',c5);

ylabel(' North-South [ km ]');  xlabel(' East-West [ km ]'); 
xlim([34.5 43]); ylim([-30 -21]);
set(gca,'xtick',15:1:55,'ytick',-38:1:-8,'tickdir','in','LineWidth',1.1,'FontSize',10,'FontName','Avenir'); 
grid on;

% print figure
figurename=sprintf('../../paperfigs/suppfig2.3_posstack_2022_300s_260129.png');
print(gcf,'-dpng','-r300',figurename);  