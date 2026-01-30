%% Stevens et al. [2015] Supplementary Info FIGURE 6 --> 2023 displacements
% GPS vs NIF: all stations data in 3 long columns
% 11 June 2021: change for flowline stations in 2012 drainage, remove NIF
% outputs
% 13 Oct 2021: remove 3\sigma outliers for the rough day boundary
% 2023 July 17: LAS revisions for reviews
% 2023 Oct 10: LAS change for 2022 lake6 array; TRACK positions
% LAS 2025-04-30: LAS add in discharge estimates for L2A, L2B from RBR :)
% 26-01-30 LAS: repository 

clear all; close all

%% load 2023 meltseason positions
load A2023_unfixedBF2_atm_datarejected_260112.mat % initial positions
station_names = {A2023.name}; sites = upper(station_names);
station_names = upper(station_names);
num_sta = length(station_names)-1; % number of stations
station_names_flip = [{'QIET'};'SQ11';'SQ12';'SQ13';'SQ14';'SQ15';'SQ16';...
    'SQ21';'SQ22';'SQ23';'SQ24';'SQ25';'SQ26';'SQ27';...
    'SQ31';'SQ32';'SQ33';'SQ34';'SQ35';'SQ36';'SQ37';'MHIH';'MLOW'];
save station_names_2023.mat station_names
open_time_L1A_2022 = 200.6;
load B2023R_BF2_UP4_sZEROm_OUT6hr_clean_260113.mat % cleaned up positions

%% positions into along- and across-flow directions based on springtime vels
springtime = [150, 160]; 
interptime_15sec = 1.736111100001381e-04; % 15-sec in decimal day
ti22_15sec = (150:interptime_15sec:160)'; % 15-second time vector

% flow angle of springtime velocity 
for i=1:1:num_sta 
    % time indices
    idx_time(i,1) = find((springtime(1)-B2023(i).neu_iu(:,1))<=0,1); % start
    idx_time(i,2) = find((springtime(2)-B2023(i).neu_iu(:,1))<=0,1); % stop

    % calculate velocity over DOY 150-160
    % horizontal vel: quick distance moved method
    bg_moved(i,1) = sqrt(((nanmean(B2023(i).neu_iu(idx_time(i,2),2))-nanmean(B2023(i).neu_iu(idx_time(i,1),2))).^2) +...
        ((nanmean(B2023(i).neu_iu(idx_time(i,2),3))-nanmean(B2023(i).neu_iu(idx_time(i,1),3))).^2)); % [ m ]
    % velocity 
    bg_moved(i,2) = bg_moved(i,1)/(springtime(2)-springtime(1)); % [ m/day ]
    bg_moved(i,3) = (bg_moved(i,1)/(springtime(2)-springtime(1)))*365.25; % [ m/yr ]
    % flow direction
    bg_moved(i,4) = (atan2d((B2023(i).neu_iu(idx_time(i,2),2)-B2023(i).neu_iu(idx_time(i,1),2)),...
        (B2023(i).neu_iu(idx_time(i,2),3)-B2023(i).neu_iu(idx_time(i,1),3)))) ;  % 0 degrees = east !!!

    % horizontal vel: linear fit to DOY 150-160 method
    B2023(i).bg_displ(:,1) = sqrt(((B2023(i).neu_iu(idx_time(i,1):idx_time(i,2),2)-B2023(i).neu_iu(idx_time(i,1),2)).^2) + ...
        ((B2023(i).neu_iu(idx_time(i,1):idx_time(i,2),3)-B2023(i).neu_iu(idx_time(i,1),3)).^2)); % [ m ]
    lin_fit(i,:) = polyfit(B2023(i).neu_iu(idx_time(i,1):idx_time(i,2),1)-ti22_15sec(1),...
        B2023(i).bg_displ(:,1), 1); 
    % velocity/flow direction -- linear fit method
    bg_moved(i,5) = B2023(i).bg_displ(end,1); % [ m ]
    bg_moved(i,6) = lin_fit(i,1); % [ m/day ]
    bg_moved(i,7) = lin_fit(i,1)*365.25; % [ m/yr ]
    if bg_moved(i,4) >= 0
        bg_moved(i,8) = 360 - bg_moved(i,4) + 90 ; % 0 degrees = north !!!
    else if bg_moved(i,4) < 0 
        bg_moved(i,8) = (-1.*bg_moved(i,4)) + 90 ; % 0 degrees = north !!!
    end
    end

    % vertical vel: linear fit to DOY 150-160 method
    B2023(i).bg_displ(:,2) = B2023(i).neu_iu(idx_time(i,1):idx_time(i,2),4)-B2023(i).neu_iu(idx_time(i,1),4); % [ m ]
    lin_fit_vertical(i,:) = polyfit(B2023(i).neu_iu(idx_time(i,1):idx_time(i,2),1)-ti22_15sec(1), B2023(i).bg_displ(:,2), 1); 
    % velocity/flow direction -- linear fit method
    bg_moved(i,9) = B2023(i).bg_displ(end,2); % [ m ]
    bg_moved(i,10) = lin_fit_vertical(i,1); % [ m/day ]
    bg_moved(i,11) = lin_fit_vertical(i,1)*365.25; % [ m/yr ]
end

%% horizontal positions into along- and across-flow directions
for i=1:1:num_sta 
   % flow angle 
   flow_angle(i,1) =  bg_moved(i,8); % flow direction [ deg ]  [0=north]
   flow_angle(i,2) =  (360-bg_moved(i,8)) + 90; % flow direction [ deg ]  [0=east]
   % displacement from t_0
   B2023(i).neu_displ(:,1) = B2023(i).neu_iu(:,1);
   B2023(i).neu_displ(:,2) = B2023(i).neu_iu(:,2) - B2023(i).neu_iu(1,2); % N
   B2023(i).neu_displ(:,3) = B2023(i).neu_iu(:,3) - B2023(i).neu_iu(1,3); % E
   B2023(i).neu_displ(:,4) = B2023(i).neu_iu(:,4) - B2023(i).neu_iu(1,4); % U
   % displacement from t_0 w/flowangle
   B2023(i).neu_displ_flow(:,1) = B2023(i).neu_iu(:,1);
   B2023(i).neu_displ_flow(:,2) = -1.*((B2023(i).neu_iu(:,2) - B2023(i).neu_iu(1,2)).*cosd(flow_angle(i,2))) + ...
       ((B2023(i).neu_iu(:,3) - B2023(i).neu_iu(1,3)).*sind(flow_angle(i,2))); % ACROSS
   B2023(i).neu_displ_flow(:,3) = ((B2023(i).neu_iu(:,3) - B2023(i).neu_iu(1,3)).*cosd(flow_angle(i,2))) + ...
       ((B2023(i).neu_iu(:,2) - B2023(i).neu_iu(1,2)).*sind(flow_angle(i,2))); % ALONG
   B2023(i).neu_displ_flow(:,4) = B2023(i).neu_iu(:,4) - B2023(i).neu_iu(1,4); % UP  
end
% check fig
% figure; clf;
% plot(B2023(14).neu_displ_flow(:,1),B2023(14).neu_displ_flow(:,2),'.'); 

%% load north lake geographic files and morlighem bed relative to North Lake
origin = [68.72, -49.53]; % M1 moulin

% load stations in polarstereographic coords
load polarstereo_stations_2023_short.mat
% polarstereo conversion needed values
radius=6378137.0;    eccen=0.08181919;    lat_true=70;    lon_posy=-45;
[moulin_x, moulin_y]=polarstereo_fwd(origin(1,1), origin(1,2),radius,eccen,lat_true,lon_posy);

% load bedmap
load('../strain_rates_2022/BMv5_for_nevis_catchment.mat'); % origin = M1 moulin 
load('../strain_rates_2022/cmapland2.mat'); % bed topo colormap  

% load lake centre points, boundaries from FASTER for 2023
% extra lake boundaries for lakes that were deemed parts of other lakes:
% load('environs_lakes_2023B_250416.mat') %  extra boundaries
% environs_lakes_2023B_extra = environs_lakes; 
%  boundaries and correct dates:
load('environs_lakes_2023B_S1S2MO_250416.mat') %  boundaries and correct dates
environs_lakes_2023B = environs_lakes; 

% load lake boundaries, centre points, drainage dates from FASTER for 2022
load('../strain_rates_2022/environs_lakes_2022B_250416.mat') % drainage dates & boundaries
environs_lakes_2022 = environs_lakes; 

% ROI that's relevant to the plotted region:
ROI_relevant_xv = [-21.75 -21.75 67.5 67.5]; % km in nevis-model space
ROI_relevant_yv = [-44 7 7 -44];
for i=1:1:length(environs_lakes_2023B.H)
    if inpolygon(environs_lakes_2023B.X_km(i),environs_lakes_2023B.Y_km(i),...
            ROI_relevant_xv,ROI_relevant_yv) == 0
        % throw a NaN
        environs_lakes_2023B.H(i) = NaN;
        environs_lakes_2023B.S(i) = NaN;
        environs_lakes_2023B.X_km(i) = NaN;
        environs_lakes_2023B.Y_km(i) = NaN;
        environs_lakes_2023B.lat(i) = NaN;
        environs_lakes_2023B.lon(i) = NaN;
        environs_lakes_2023B.drainage_type_num(i) = NaN;
        environs_lakes_2023B.laketypeing_dates(i,1:7) = [NaN, NaN, NaN, NaN, NaN, NaN, NaN];
    else
    end
end
% redefine lake_typeing to just include lakes within nevis domain
lake_typeing = environs_lakes_2023B.laketypeing_dates; % after throwing out lakes on the nevis domain border
% polygon ROI check
figure(2); plot(ROI_relevant_xv,ROI_relevant_yv)
hold on; plot(environs_lakes_2023B.X_km,environs_lakes_2023B.Y_km,'o');
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

%% load pressure logger records
RB32_L2A_pos = [68.62908,-49.00388];
[RB32_L2A_x,RB32_L2A_y]=polarstereo_fwd(68.62908,-49.00388,radius,eccen,lat_true,lon_posy); % [m]
RB32_L2A_x_moulin = (RB32_L2A_x-moulin_x)./1e3;
RB32_L2A_y_moulin = (RB32_L2A_y-moulin_y)./1e3;

RB42_L2B_pos = [68.57236,-49.08983];
[RB42_L2B_x,RB42_L2B_y]=polarstereo_fwd(68.57236,-49.08983,radius,eccen,lat_true,lon_posy); % [m]
RB42_L2B_x_moulin = (RB42_L2B_x-moulin_x)./1e3;
RB42_L2B_y_moulin = (RB42_L2B_y-moulin_y)./1e3;

% cleaned up and discharge calculated:
load('../strain_rates_2022/volume/RB32_2023_interp_volume_discharge.mat')
load('../strain_rates_2022/volume/RB42_2023_interp_volume_discharge.mat')

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

%% Figure 200s - C6 STATIC
figure(2); clf;
set(gcf,'Units','centimeters','Position',1.5.*[0.5 0 19.3 20]);
    
xmin = 202.4; xmax = 203.3; 
ymin = -0.25; ymax = 4.5;
offset_all = 7.5;

    axeRBR = axes('Position',[0.059 0.80 0.29 0.1725],'Box','on','NextPlot','add');    
    axe0 = axes('Position',[0.07+0.29+0.08 0.79 0.195 0.1825],'Box','on','NextPlot','add');

    axe1 = axes('Position',[0.059 0.04 0.29 0.94],'Box','on','NextPlot','add');
    ylabel('GNSS Station','FontSize',11);
    set(gca, 'XTickLabelMode', 'auto');
    set(gca,'YTick',0:0.25:ymax,'YTickLabel',{'SQ27';'';'SQ26';'';'SQ25';'';'SQ24';'';'SQ23';'';'SQ22';...
        '';'SQ21';'';'';'';'';'';'';'';''})
    
    axe2 = axes('Position',[0.055+0.30 0.04 0.29 0.94],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'auto');
    set(gca,'YTick',0:0.25:ymax,'YTickLabel',[]);
    xlabel('Day of Year, 2023  [ UTC ]','FontSize',11,'Fontname','Avenir')

    axe3 = axes('Position',[0.053+0.30+0.30 0.04 0.29 0.94],'Box','on','NextPlot','add');
    ylabel('Displacement  [ m ]','FontSize',11,'Fontname','Avenir');
    set(gca, 'XTickLabelMode', 'auto');
    set(gca,'yaxislocation','right','ytick',0:0.25:ymax)
    
    TriangleSize = 7;
    offset = 0.5; % displacements vertical offset
    offset_v = offset.*(0:1:21); % displacements y-intercept

% begin GPS displacements
for i=10:16 % 1150s num_sta

    % find within the halfhour of xmin (over a 1-hr window)
    time_before = 1/24; % [decimal days]
    [II] = find(B2023(i).neu_displ_flow(2:end,1)>(xmin-(time_before/2)) & ...
          B2023(i).neu_displ_flow(2:end,1)<(xmin+(time_before/2)));    
    pos_30min_before_open_time_L1A_2022(i,:) = nanmean(B2023(i).neu_displ_flow(II,1:4));
    % find background velocity 1 day before xmin (over a 24-hr window)
    time_before = 0.2; % [decimal days]
    [I] = find(B2023(i).neu_displ_flow(2:end,1)>(xmin-(time_before)) & ...
          B2023(i).neu_displ_flow(2:end,1)<(xmin+0.1));
    vel_30min_before(i,1) = B2023(i).neu_displ_flow(I(1),1); % time
        % slopes
        s_across = polyfit(B2023(i).neu_displ_flow(I,1),B2023(i).neu_displ_flow(I,2),1); % slope in across [m/day]
        s_along = polyfit(B2023(i).neu_displ_flow(I,1),B2023(i).neu_displ_flow(I,3),1); % slope in along [m/day]
        s_up = polyfit(B2023(i).neu_displ_flow(I,1),B2023(i).neu_displ_flow(I,4),1); % slope in up [m/day]
        % velocities
        vel_30min_before(i,2) = s_along(1); % [ m/d ]
        vel_30min_before(i,3) = s_across(1); % [ m/d ]
        vel_30min_before(i,4) = s_up(1); % [ m/d ]

  axes(axe1);% E
  title('Flowline (263-270˚) displacement [ m ]','Fontname','Avenir','FontSize',11);
  plot_data_along = offset_all + (B2023(i).neu_displ_flow(:,3)-pos_30min_before_open_time_L1A_2022(i,3)) - ...
                    offset_v(i) - (vel_30min_before(i,2).*(B2023(i).neu_displ_flow(:,1)-xmin));
  plot(B2023(i).neu_displ_flow(:,1),plot_data_along,'.','markerSize',2.5); 
  hold on; grid on; set(gca,'XMinorGrid','on','YMinorGrid','on','LineWidth',1.1,'FontSize',10,'FontName','Avenir')
  text(xmin+0.03, 3.4, 'a','FontWeight','bold','FontSize',12);
  xlim([xmin xmax]); ylim([ymin ymax])
  box on

  axes(axe2) %N
  title('Across-flow (353-360˚) displacement [ m ]','Fontname','Avenir','FontSize',11)
  plot_data_across = offset_all + (B2023(i).neu_displ_flow(:,2)-pos_30min_before_open_time_L1A_2022(i,2)) - ...
      offset_v(i) - (vel_30min_before(i,3).*(B2023(i).neu_displ_flow(:,1)-xmin));
  plot(B2023(i).neu_displ_flow(:,1),plot_data_across,'.','markerSize',2.5); 
  hold on; grid on; set(gca,'XMinorGrid','on','YMinorGrid','on','LineWidth',1.1,'FontSize',10,'FontName','Avenir')
  text(xmin+0.03, ymax-0.07, 'b','FontWeight','bold','FontSize',12);
  xlim([xmin xmax]); ylim([ymin ymax])
  box on
  
  axes(axe3) %U
  title('Vertical displacement [ m ]','Fontname','Avenir','FontSize',11)
  plot_data_up = offset_all + (B2023(i).neu_displ_flow(:,4)-pos_30min_before_open_time_L1A_2022(i,4)) - ...
      offset_v(i) - (vel_30min_before(i,4).*(B2023(i).neu_displ_flow(:,1)-xmin));
  plot(B2023(i).neu_displ_flow(:,1),plot_data_up,'.','markerSize',2.5); 
  hold on; grid on; set(gca,'XMinorGrid','on','YMinorGrid','on','LineWidth',1.1,'FontSize',10,'FontName','Avenir')
  text(xmin+0.03, ymax-0.07, 'c','FontWeight','bold','FontSize',12);
  xlim([xmin xmax]); ylim([ymin ymax])

end

% % DISPLACEMENTS for just MLOW
% for i=2
%     % find within the halfhour of xmin (over a 1-hr window)
%     open_time_C6_2023 = 200.6;
%     time_before = 1/24; % [decimal days]
%     [II] = find(B2023(i).neu_displ_flow(2:end,1)>(xmin-(time_before/2)) & ...
%           B2023(i).neu_displ_flow(2:end,1)<(xmin+(time_before/2)));    
%     pos_30min_before_open_time_L1A_2022(i,:) = nanmean(B2023(i).neu_displ_flow(II,1:4));
%     % find background velocity 1 day before xmin (over a 24-hr window)
%     time_before = 0.25; % [decimal days]
%     [I] = find(B2023(i).neu_displ_flow(2:end,1)>(xmin-(time_before)) & ...
%           B2023(i).neu_displ_flow(2:end,1)<(xmin+0.1));
%     vel_30min_before(i,1) = B2023(i).neu_displ_flow(I(1),1); % time
%         % slopes
%         s_across = polyfit(B2023(i).neu_displ_flow(I,1),B2023(i).neu_displ_flow(I,2),1); % slope in across [m/day]
%         s_along = polyfit(B2023(i).neu_displ_flow(I,1),B2023(i).neu_displ_flow(I,3),1); % slope in along [m/day]
%         s_up = polyfit(B2023(i).neu_displ_flow(I,1),B2023(i).neu_displ_flow(I,4),1); % slope in up [m/day]
%         % velocities
%         vel_30min_before(i,2) = s_along(1); % [ m/d ]
%         vel_30min_before(i,3) = s_across(1); % [ m/d ]
%         vel_30min_before(i,4) = s_up(1); % [ m/d ]
% 
%   axes(axe1);% E
%   title('Flowline (277-279?) displacement [ m ]','Fontname','Avenir','FontSize',11);
%   plot_data_along = offset_all + (B2023(i).neu_displ_flow(:,3)-pos_30min_before_open_time_L1A_2022(i,3)) - ...
%                     offset_v(i+7) - (vel_30min_before(i,2).*(B2023(i).neu_displ_flow(:,1)-xmin));
%   plot(B2023(i).neu_displ_flow(:,1),plot_data_along,'.','markerSize',2.5); 
%   hold on; grid on; set(gca,'XMinorGrid','on','YMinorGrid','on','LineWidth',1.1,'FontSize',10,'FontName','Avenir')
%   text(xmin+0.03, ymax-0.07, 'a','FontWeight','bold','FontSize',12);
%   xlim([xmin xmax]); ylim([ymin ymax])
%   box on
% 
%   axes(axe2) %N
%   title('Across-flow (7-9?) displacement [ m ]','Fontname','Avenir','FontSize',11)
%   plot_data_across = offset_all + (B2023(i).neu_displ_flow(:,2)-pos_30min_before_open_time_L1A_2022(i,2)) - ...
%       offset_v(i+7) - (vel_30min_before(i,3).*(B2023(i).neu_displ_flow(:,1)-xmin));
%   plot(B2023(i).neu_displ_flow(:,1),plot_data_across,'.','markerSize',2.5); 
%   hold on; grid on; set(gca,'XMinorGrid','on','YMinorGrid','on','LineWidth',1.1,'FontSize',10,'FontName','Avenir')
%   text(xmin+0.03, ymax-0.07, 'b','FontWeight','bold','FontSize',12);
%   xlim([xmin xmax]); ylim([ymin ymax])
%   box on
% 
%   axes(axe3) %U
%   title('Vertical displacement [ m ]','Fontname','Avenir','FontSize',11)
%   plot_data_up = offset_all + (B2023(i).neu_displ_flow(:,4)-pos_30min_before_open_time_L1A_2022(i,4)) - ...
%       offset_v(i+7) - (vel_30min_before(i,4).*(B2023(i).neu_displ_flow(:,1)-xmin));
%   plot(B2023(i).neu_displ_flow(:,1),plot_data_up,'.','markerSize',2.5); 
%   hold on; grid on; set(gca,'XMinorGrid','on','YMinorGrid','on','LineWidth',1.1,'FontSize',10,'FontName','Avenir')
%   text(xmin+0.03, ymax-0.07, 'c','FontWeight','bold','FontSize',12);
%   xlim([xmin xmax]); ylim([ymin ymax])
% 
% end


% RBR
axes(axeRBR)
circle_size = 2;
text(xmin+0.02, -350, 'd','FontWeight','bold','FontSize',11); hold on;
%if i==10; rectangle('Position',[time_real -1 time_vector(k+1)-time_real 12],'FaceColor',[0.8 0.8 0.8],'EdgeColor','none'); end
plot(RB32_2023_interp_volume_discharge.time_dDAY,...
    RB32_2023_interp_volume_discharge.discharge_dDAY_RB32_smooth,...
    '-','Color',[1.0 0.6 0.6],'LineWidth',1.2); hold on
plot(RB42_2023_interp_volume_discharge.time_dDAY,...
    RB42_2023_interp_volume_discharge.discharge_dDAY_RB42_smooth,...
    '-','Color',[0.6 0.6 1.0],'LineWidth',1.2); hold on
%plot(RB32_2023(:,1), RB32_2023(:,3)+1.29476,'.','Color',[1.0 0.6 0.6],'MarkerSize',circle_size*2); hold on;
%plot(RB42_2023(:,1), RB42_2023(:,3)+0.973573,'.','Color',[0.6 0.6 1.0],'MarkerSize',circle_size*2);
xlim([xmin xmax]); ylim([-5000 100]); %ylim([-1 11]);
text(202.7, -3500, {'L2A';'(RB32)'},'FontName','AvenirBold','FontSize',14,'FontAngle','italic','Color',[1.0 0.6 0.6],'HorizontalAlignment','center')
text(203.0, -1500, {'L2B';'(RB42)'},'FontName','AvenirBold','FontSize',14,'FontAngle','italic','Color',[0.6 0.6 1.0],'HorizontalAlignment','center')
ylabel('Lake Discharge [ m^{3} s^{-1} ]','FontSize',10,'FontName','Avenir');
hold on; grid on; 
set(gca,'XMinorGrid','on','YMinorGrid','on','LineWidth',1.1,'FontSize',10,...
    'FontName','Avenir','XTickLabel',[])

    
axes(axe0);
hold on;
text(13.25, -10.5, 'e','FontWeight','bold','FontSize',11);
hold on;
% plot(xy_sta_23(:,1),xy_sta_23(:,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');
% RBRs
scatter(RB32_L2A_x(:,1)./1e3,RB32_L2A_y(:,1)./1e3,78,'bp','filled','MarkerEdgeColor','none','MarkerFaceColor',[1.0 0.6 0.6])
scatter(RB42_L2B_x(:,1)./1e3,RB42_L2B_y(:,1)./1e3,78,'bp','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 1.0])% hydrofractures
% lake names
text(20.15, -11.25, 'L2A','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(19.5, -15.5, 'L2B','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(RB32_L2A_x(:,1)./1e3+0.25,RB32_L2A_y(:,1)./1e3,'RB32','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(RB42_L2B_x(:,1)./1e3+0.25,RB42_L2B_y(:,1)./1e3,'RB42','FontSize',11,'FontName','Avenir','FontAngle','italic')

% lakes: FASTER S2 2023
[HF_lakes_to_plot] = find(environs_lakes_2023B.laketypeing_dates(:,7)==1); % just HF lakes  
[MU_lakes_to_plot] = find(environs_lakes_2023B.laketypeing_dates(:,7)==2); % just moulin-draining lakes 
[lakes_to_plot] = find(~isnan(environs_lakes_2023B.laketypeing_dates(:,7))); % all good lakes   
[NE_lakes_to_plot] = find(environs_lakes_2023B.laketypeing_dates(:,7)==5); % just no-exit freezers
for i=1:1:length(lakes_to_plot)
    fill3(environs_lakes_2023B.boundaries(lakes_to_plot(i)).XY_km_local(:,1), ...
        environs_lakes_2023B.boundaries(lakes_to_plot(i)).XY_km_local(:,2), ...
        -3000.*ones(length(environs_lakes_2023B.boundaries(lakes_to_plot(i)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',lake_blue(1:3));
end
for i=1:1:length(HF_lakes_to_plot)
    fill3(environs_lakes_2023B.boundaries(HF_lakes_to_plot(i)).XY_km_local(:,1), ...
        environs_lakes_2023B.boundaries(HF_lakes_to_plot(i)).XY_km_local(:,2), ...
        -3000.*ones(length(environs_lakes_2023B.boundaries(HF_lakes_to_plot(i)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',goldenrod);
end
for i=1:1:length(MU_lakes_to_plot)
    fill3(environs_lakes_2023B.boundaries(MU_lakes_to_plot(i)).XY_km_local(:,1), ...
        environs_lakes_2023B.boundaries(MU_lakes_to_plot(i)).XY_km_local(:,2), ...
        -3000.*ones(length(environs_lakes_2023B.boundaries(MU_lakes_to_plot(i)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',[0.4 0.8 0.8]);
end
for i=1:1:length(NE_lakes_to_plot) % No-exit freezers 
    fill3(environs_lakes_2023B.boundaries(NE_lakes_to_plot(i)).XY_km_local(:,1), ...
        environs_lakes_2023B.boundaries(NE_lakes_to_plot(i)).XY_km_local(:,2), ...
        -3000.*ones(length(environs_lakes_2023B.boundaries(NE_lakes_to_plot(i)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',plum(1:3));
end
% % L1A from 2022 outline (2022 ID = 81)
% fill3(environs_lakes_2022.boundaries(81).XY_km_local(:,1), ...
%         environs_lakes_2022.boundaries(81).XY_km_local(:,2), ...
%         -3000.*ones(length(environs_lakes_2022.boundaries(81).XY_km_local(:,2)),1),...
%          'o','MarkerSize',1,'EdgeColor','none','FaceColor',goldenrod);
% % L1C from 2022 outline (2022 ID = 77)
% fill3(environs_lakes_2022.boundaries(77).XY_km_local(:,1), ...
%         environs_lakes_2022.boundaries(77).XY_km_local(:,2), ...
%         -3000.*ones(length(environs_lakes_2022.boundaries(77).XY_km_local(:,2)),1),...
%          'o','MarkerSize',1,'EdgeColor','none','FaceColor',goldenrod);
% % L2A from 2022 outline (2022 ID = 156)
% fill3(environs_lakes_2022.boundaries(156).XY_km_local(:,1), ...
%         environs_lakes_2022.boundaries(156).XY_km_local(:,2), ...
%         -3000.*ones(length(environs_lakes_2022.boundaries(156).XY_km_local(:,2)),1),...
%          'o','MarkerSize',1,'EdgeColor','none','FaceColor',goldenrod);

% GNSS stations
text(xy_sta_23_short(10:16,1)-1.20,xy_sta_23_short(10:16,2)-0.5,station_names(10:16),'FontSize',9,'Rotation',0);
plot(xy_sta_23_short(10,1),xy_sta_23_short(10,2),'k^','MarkerSize',TriangleSize-0,'MarkerFaceColor',c1);
plot(xy_sta_23_short(11,1),xy_sta_23_short(11,2),'k^','MarkerSize',TriangleSize-0,'MarkerFaceColor',c2);
plot(xy_sta_23_short(12,1),xy_sta_23_short(12,2),'k^','MarkerSize',TriangleSize-0,'MarkerFaceColor',c3);
plot(xy_sta_23_short(13,1),xy_sta_23_short(13,2),'k^','MarkerSize',TriangleSize-0,'MarkerFaceColor',c4);
plot(xy_sta_23_short(14,1),xy_sta_23_short(14,2),'k^','MarkerSize',TriangleSize-0,'MarkerFaceColor',c5);
plot(xy_sta_23_short(15,1),xy_sta_23_short(15,2),'k^','MarkerSize',TriangleSize-0,'MarkerFaceColor',c6);
plot(xy_sta_23_short(16,1),xy_sta_23_short(16,2),'k^','MarkerSize',TriangleSize-0,'MarkerFaceColor',c7);

ylabel(' North-South [ km ]');  xlabel(' East-West [ km ]'); 
xlim([13 23]); ylim([-21 -10]);
set(gca,'xtick',10:1:55,'ytick',-38:1:-8,'tickdir','in','LineWidth',1.1,'FontSize',8,'FontName','Avenir'); 
grid on;

% print figure
% figurename=sprintf('type_cases/C5/pos_stack_sZEROm_200s_250430.png');
% print(gcf,'-dpng','-r600',figurename);  

%% Figure 2023 200s - MOVIE
colors = colororder;
time_vector = xmin:(1/48):xmax; % half-hour windows

for k = 23
% for k = 23:length(time_vector)-1
    time_real = time_vector(k); % real time

% Figure 200s MOVIE
figure(3); clf;
set(gcf,'Units','centimeters','Position',1.5.*[0.5 0 19.3 20]);
    
xmin = 202.4; xmax = 203.3; 
ymin = -0.25; ymax = 4.5;
offset_all = 7.5;
    
    axeRBR = axes('Position',[0.059 0.80 0.29 0.1725],'Box','on','NextPlot','add'); 
    axe0 = axes('Position',[0.07+0.29+0.08 0.79 0.195 0.1825],'Box','on','NextPlot','add');
    axe00 = axes('Position',[0.07+0.29+0.29+0.08 0.755 0.200 0.2025],'Box','on',...
        'NextPlot','add','YTickLabel',[],'XTickLabel',[]);

    axe1 = axes('Position',[0.059 0.04 0.29 0.94],'Box','on','NextPlot','add');
    ylabel('GNSS Station','FontSize',11);
    set(gca, 'XTickLabelMode', 'auto');
    set(gca,'YTick',0:0.25:ymax,'YTickLabel',{'SQ27';'';'SQ26';'';'SQ25';'';'SQ24';'';'SQ23';'';'SQ22';...
        '';'SQ21';'';'';'';'';'';'';'';''})
    
    axe2 = axes('Position',[0.055+0.30 0.04 0.29 0.94],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'auto');
    set(gca,'YTick',0:0.25:ymax,'YTickLabel',[]);
    xlabel('Day of Year, 2023  [ UTC ]','FontSize',11,'Fontname','Avenir')

    axe3 = axes('Position',[0.053+0.30+0.30 0.04 0.29 0.94],'Box','on','NextPlot','add');
    ylabel('Displacement  [ m ]','FontSize',11,'Fontname','Avenir');
    set(gca, 'XTickLabelMode', 'auto');
    set(gca,'yaxislocation','right','ytick',0:0.25:ymax)
    
    TriangleSize = 7;
    offset = 0.5; % displacements vertical offset
    offset_v = offset.*(0:1:21); % displacements y-intercept

% begin GPS displacements
for i=10:16 % 1150s num_sta
    % find within the halfhour of xmin (over a 1-hr window)
    open_time_C6_2023 = 202.7;
    time_before = 1/24; % [decimal days]
    [II] = find(B2023(i).neu_displ_flow(2:end,1)>(xmin-(time_before/2)) & ...
          B2023(i).neu_displ_flow(2:end,1)<(xmin+(time_before/2)));    
    pos_30min_before_open_time_L1A_2022(i,:) = nanmean(B2023(i).neu_displ_flow(II,1:4));
    % find background velocity 1 day before xmin (over a 24-hr window)
    time_before = 0.2; % [decimal days]
    [I] = find(B2023(i).neu_displ_flow(2:end,1)>(xmin-(time_before)) & ...
          B2023(i).neu_displ_flow(2:end,1)<(xmin+0.1));
    vel_30min_before(i,1) = B2023(i).neu_displ_flow(I(1),1); % time
        % slopes
        s_across = polyfit(B2023(i).neu_displ_flow(I,1),B2023(i).neu_displ_flow(I,2),1); % slope in across [m/day]
        s_along = polyfit(B2023(i).neu_displ_flow(I,1),B2023(i).neu_displ_flow(I,3),1); % slope in along [m/day]
        s_up = polyfit(B2023(i).neu_displ_flow(I,1),B2023(i).neu_displ_flow(I,4),1); % slope in up [m/day]
        % velocities
        vel_30min_before(i,2) = s_along(1); % [ m/d ]
        vel_30min_before(i,3) = s_across(1); % [ m/d ]
        vel_30min_before(i,4) = s_up(1); % [ m/d ]

  axes(axe1);% E
  title('Flowline (263-270˚) displacement [ m ]','Fontname','Avenir','FontSize',11);
  hold on; 
  if i==10; rectangle('Position',[time_real -100 time_vector(k+1)-time_real 200],'FaceColor',[0.8 0.8 0.8],'EdgeColor','none'); end
  plot_data_along = offset_all + (B2023(i).neu_displ_flow(:,3)-pos_30min_before_open_time_L1A_2022(i,3)) - ...
                    offset_v(i) - (vel_30min_before(i,2).*(B2023(i).neu_displ_flow(:,1)-xmin));
  plot(B2023(i).neu_displ_flow(:,1),plot_data_along,'.','markerSize',2.5,'color',colors(i-9,:)); 
  grid on; set(gca,'XMinorGrid','on','YMinorGrid','on','LineWidth',1.1,'FontSize',10,'FontName','Avenir')
  text(xmin+0.03, 3.4, 'a','FontWeight','bold','FontSize',12);
  xlim([xmin xmax]); ylim([ymin ymax])
  box on

  axes(axe2) %N
  title('Across-flow (353-360˚) displacement [ m ]','Fontname','Avenir','FontSize',11)
  hold on; 
  if i==10 
      rectangle('Position',[time_real -100 time_vector(k+1)-time_real 200],'FaceColor',[0.8 0.8 0.8],'EdgeColor','none'); 
  end
  plot_data_across = offset_all + (B2023(i).neu_displ_flow(:,2)-pos_30min_before_open_time_L1A_2022(i,2)) - ...
      offset_v(i) - (vel_30min_before(i,3).*(B2023(i).neu_displ_flow(:,1)-xmin));
  plot(B2023(i).neu_displ_flow(:,1),plot_data_across,'.','markerSize',2.5,'color',colors(i-9,:)); 
  hold on; grid on; set(gca,'XMinorGrid','on','YMinorGrid','on','LineWidth',1.1,'FontSize',10,'FontName','Avenir')
  text(xmin+0.03, ymax-0.07, 'b','FontWeight','bold','FontSize',12);
  xlim([xmin xmax]); ylim([ymin ymax])
  box on
  
  axes(axe3) %U
  title('Vertical displacement [ m ]','Fontname','Avenir','FontSize',11)
  hold on; 
  if i==10 
      rectangle('Position',[time_real -100 time_vector(k+1)-time_real 200],'FaceColor',[0.8 0.8 0.8],'EdgeColor','none'); 
      rectangle('Position',[xmin+0.01 3.35 xmax-xmin-0.02 1.14],'FaceColor',[1 1 1],'EdgeColor','none'); 
  end
  plot_data_up = offset_all + (B2023(i).neu_displ_flow(:,4)-pos_30min_before_open_time_L1A_2022(i,4)) - ...
      offset_v(i) - (vel_30min_before(i,4).*(B2023(i).neu_displ_flow(:,1)-xmin));
  plot(B2023(i).neu_displ_flow(:,1),plot_data_up,'.','markerSize',2.5,'color',colors(i-9,:)); 
  hold on; grid on; set(gca,'XMinorGrid','on','YMinorGrid','on','LineWidth',1.1,'FontSize',10,'FontName','Avenir')
  text(xmin+0.03, 3.20, 'c','FontWeight','bold','FontSize',12);
  xlim([xmin xmax]); ylim([ymin ymax])
  
  % calculate east for each station
  pos_30min_before_east(i,:) = nanmean(B2023(i).neu_displ(II,1:4));
  vel_30min_before_east(i,:) = B2023(i).neu_displ(I(1),1); % time
        % slopes
        s_across_east = polyfit(B2023(i).neu_displ(I,1),B2023(i).neu_displ(I,2),1); % slope in across [m/day]
        s_along_east = polyfit(B2023(i).neu_displ(I,1),B2023(i).neu_displ(I,3),1); % slope in along [m/day]
        s_up_east = polyfit(B2023(i).neu_displ(I,1),B2023(i).neu_displ(I,4),1); % slope in up [m/day]
        % velocities
        vel_30min_before_east(i,2) = s_along_east(1); % [ m/d ] E
        vel_30min_before_east(i,3) = s_across_east(1); % [ m/d ] N
        vel_30min_before_east(i,4) = s_up_east(1); % [ m/d ] U
  % saving horizontals for each station
  B2023(i).east_displacement_time = (B2023(i).neu_displ(:,3)-pos_30min_before_east(i,3)) - ...
                   (vel_30min_before_east(i,2).*(B2023(i).neu_displ(:,1)-xmin));
  B2023(i).north_displacement_time = (B2023(i).neu_displ(:,2)-pos_30min_before_east(i,2)) - ...
                   (vel_30min_before_east(i,3).*(B2023(i).neu_displ(:,1)-xmin));
  % saving vertical displacement for each station
  B2023(i).vert_displacement_time = (B2023(i).neu_displ(:,4)-pos_30min_before_east(i,4)) - ...
      (vel_30min_before_east(i,4).*(B2023(i).neu_displ(:,1)-xmin));
end 

% RBR
axes(axeRBR)
circle_size = 2;
text(xmin+0.02, -350, 'd','FontWeight','bold','FontSize',11); hold on;
rectangle('Position',[time_real -5000 time_vector(k+1)-time_real 5099],'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
plot(RB32_2023_interp_volume_discharge.time_dDAY,...
    RB32_2023_interp_volume_discharge.discharge_dDAY_RB32_smooth,...
    '-','Color','k','LineWidth',1.2); hold on
plot(RB42_2023_interp_volume_discharge.time_dDAY,...
    RB42_2023_interp_volume_discharge.discharge_dDAY_RB42_smooth,...
    '-','Color',lake_blue(1:3),'LineWidth',1.2); hold on
xlim([xmin xmax]); ylim([-5000 100]); %ylim([-1 11]);
text(202.7, -3500, {'L2A';'(RB32)'},'FontName','AvenirBold','FontSize',14,'FontAngle','italic','Color','k','HorizontalAlignment','center')
text(203.0, -1500, {'L2B';'(RB42)'},'FontName','AvenirBold','FontSize',14,'FontAngle','italic','Color',lake_blue(1:3),'HorizontalAlignment','center')
ylabel('Lake Discharge [ m^{3} s^{-1} ]','FontSize',9,'FontName','Avenir');
hold on; grid on; 
set(gca,'XMinorGrid','on','YMinorGrid','on','LineWidth',1.1,'FontSize',9,...
    'FontName','Avenir','XTickLabel',[])
  

for j=10:16 % loop over stations for this time of the movie
    % indices for time:
    [III] = find(B2023(j).neu_displ(2:end,1)>(time_vector(k)) & ...
          B2023(j).neu_displ(2:end,1)<(time_vector(k+1)));
    % vertical displacement at the right time of the movie:
    plot_vert_displacement_time(j,1) = nanmean(B2023(j).vert_displacement_time(III));
    % horiztonal displacement
    plot_north_displacement_time(j,1) = B2023(j).north_displacement_time(III(end))-B2023(j).north_displacement_time(III(1)); % north [ m ] 
    plot_east_displacement_time(j,1) = B2023(j).east_displacement_time(III(end))-B2023(j).east_displacement_time(III(1)); % east [ m ]
end
 

axes(axe00);
text(12.35, -9.5, 'e','FontWeight','bold','FontSize',11);
hold on;
% RBR
scatter(RB32_L2A_x(:,1)./1e3,RB32_L2A_y(:,1)./1e3,78,'bp','filled','MarkerEdgeColor','none','MarkerFaceColor',[1.0 0.6 0.6])
scatter(RB42_L2B_x(:,1)./1e3,RB42_L2B_y(:,1)./1e3,78,'bp','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 1.0])
hold on;
% quiver (horizontals)
quiver_boost = 17;
quiver(xy_sta_23_short(10,1),xy_sta_23_short(10,2),...
    plot_east_displacement_time(10)*quiver_boost,plot_north_displacement_time(10)*quiver_boost,...
    "LineWidth",1.3,"Color",c1,"AutoScale","off",'MaxHeadSize',0.5);
quiver(xy_sta_23_short(11,1),xy_sta_23_short(11,2),...
    plot_east_displacement_time(11)*quiver_boost,plot_north_displacement_time(11)*quiver_boost,...
    "LineWidth",1.3,"Color",c2,"AutoScale","off",'MaxHeadSize',0.5);
quiver(xy_sta_23_short(12,1),xy_sta_23_short(12,2),...
    plot_east_displacement_time(12)*quiver_boost,plot_north_displacement_time(12)*quiver_boost,...
    "LineWidth",1.3,"Color",c3,"AutoScale","off",'MaxHeadSize',0.5);
quiver(xy_sta_23_short(13,1),xy_sta_23_short(13,2),...
    plot_east_displacement_time(13)*quiver_boost,plot_north_displacement_time(13)*quiver_boost,...
    "LineWidth",1.3,"Color",c4,"AutoScale","off",'MaxHeadSize',0.5);
quiver(xy_sta_23_short(14,1),xy_sta_23_short(14,2),...
    plot_east_displacement_time(14)*quiver_boost,plot_north_displacement_time(14)*quiver_boost,...
    "LineWidth",1.3,"Color",c5,"AutoScale","off",'MaxHeadSize',0.5);
quiver(xy_sta_23_short(15,1),xy_sta_23_short(15,2),...
    plot_east_displacement_time(15)*quiver_boost,plot_north_displacement_time(15)*quiver_boost,...
    "LineWidth",1.3,"Color",c6,"AutoScale","off",'MaxHeadSize',0.5);
quiver(xy_sta_23_short(16,1),xy_sta_23_short(16,2),...
    plot_east_displacement_time(16)*quiver_boost,plot_north_displacement_time(16)*quiver_boost,...
    "LineWidth",1.3,"Color",c7,"AutoScale","off",'MaxHeadSize',0.5);
% quiver scalebar
quiver(22.25,-19.75,-0.10*quiver_boost,0*quiver_boost,...
    "LineWidth",1.3,"Color",'k',"AutoScale","off",'MaxHeadSize',0.6);
text(20.5, -19.0, '10 cm','FontSize',10,'FontName','AvenirBold')

% lakes: FASTER S2 2023
[HF_lakes_to_plot] = find(environs_lakes_2023B.laketypeing_dates(:,7)==1); % just HF lakes  
[MU_lakes_to_plot] = find(environs_lakes_2023B.laketypeing_dates(:,7)==2); % just moulin-draining lakes 
[lakes_to_plot] = find(~isnan(environs_lakes_2023B.laketypeing_dates(:,7))); % all good lakes   
[NE_lakes_to_plot] = find(environs_lakes_2023B.laketypeing_dates(:,7)==5); % just no-exit freezers
for i=1:1:length(lakes_to_plot)
    fill3(environs_lakes_2023B.boundaries(lakes_to_plot(i)).XY_km_local(:,1), ...
        environs_lakes_2023B.boundaries(lakes_to_plot(i)).XY_km_local(:,2), ...
        -3000.*ones(length(environs_lakes_2023B.boundaries(lakes_to_plot(i)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',lake_blue(1:3));
end
for i=1:1:length(HF_lakes_to_plot)
    fill3(environs_lakes_2023B.boundaries(HF_lakes_to_plot(i)).XY_km_local(:,1), ...
        environs_lakes_2023B.boundaries(HF_lakes_to_plot(i)).XY_km_local(:,2), ...
        -3000.*ones(length(environs_lakes_2023B.boundaries(HF_lakes_to_plot(i)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',goldenrod);
end
for i=1:1:length(MU_lakes_to_plot)
    fill3(environs_lakes_2023B.boundaries(MU_lakes_to_plot(i)).XY_km_local(:,1), ...
        environs_lakes_2023B.boundaries(MU_lakes_to_plot(i)).XY_km_local(:,2), ...
        -3000.*ones(length(environs_lakes_2023B.boundaries(MU_lakes_to_plot(i)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',[0.4 0.8 0.8]);
end
for i=1:1:length(NE_lakes_to_plot) % No-exit freezers 
    fill3(environs_lakes_2023B.boundaries(NE_lakes_to_plot(i)).XY_km_local(:,1), ...
        environs_lakes_2023B.boundaries(NE_lakes_to_plot(i)).XY_km_local(:,2), ...
        -3000.*ones(length(environs_lakes_2023B.boundaries(NE_lakes_to_plot(i)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',plum(1:3));
end

% GNSS stations
text(xy_sta_23_short(10:16,1)-0.5,xy_sta_23_short(10:16,2)+0.75,station_names(10:16),'FontSize',9,'Rotation',0);
% GNSS stations uplift
scatter(xy_sta_23_short(10,1),xy_sta_23_short(10,2),TriangleSize+58,plot_vert_displacement_time(10,1),'filled','^');
scatter(xy_sta_23_short(11,1),xy_sta_23_short(11,2),TriangleSize+58,plot_vert_displacement_time(11,1),'filled','^');
scatter(xy_sta_23_short(12,1),xy_sta_23_short(12,2),TriangleSize+58,plot_vert_displacement_time(12,1),'filled','^');
scatter(xy_sta_23_short(13,1),xy_sta_23_short(13,2),TriangleSize+58,plot_vert_displacement_time(13,1),'filled','^');
scatter(xy_sta_23_short(14,1),xy_sta_23_short(14,2),TriangleSize+58,plot_vert_displacement_time(14,1),'filled','^');
scatter(xy_sta_23_short(15,1),xy_sta_23_short(15,2),TriangleSize+58,plot_vert_displacement_time(15,1),'filled','^');
scatter(xy_sta_23_short(16,1),xy_sta_23_short(16,2),TriangleSize+58,plot_vert_displacement_time(16,1),'filled','^');
% RBR
scatter(RB32_L2A_x_moulin(:,1),RB32_L2A_y_moulin(:,1),98,'bp','filled','MarkerEdgeColor','none','MarkerFaceColor','k')
scatter(RB42_L2B_x_moulin(:,1),RB42_L2B_y_moulin(:,1),98,'bp','filled','MarkerEdgeColor','none','MarkerFaceColor',lake_blue(1:3))
text(RB32_L2A_x_moulin(:,1)+0.35,RB32_L2A_y_moulin(:,1),'RB32','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(RB42_L2B_x_moulin(:,1)+0.35,RB42_L2B_y_moulin(:,1)+0.25,'RB42','FontSize',11,'FontName','Avenir','FontAngle','italic')
% lake names
text(17.15, -12.5, 'L2A','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(14.0, -16.5, 'L2B','FontSize',11,'FontName','Avenir','FontAngle','italic')
% colorbar
colormap('jet')
clb = colorbar('Position',[0.71 0.772 0.007 0.18]); clim([0 0.75]);
set(clb,'FontName','Avenir','FontSize',10,'YTick',0:0.25:0.75); 
set(get(clb,'ylabel'),'String',{'Uplift  [ m ]'},'FontSize',11,'Fontname','AvenirBold');
title(sprintf('2023/%6.3f--%6.3f',time_real,time_vector(k+1)),'FontName','Avenir','FontSize',12,'FontWeight','bold')
%ylabel('North-South [ km ]');  xlabel('East-West [ km ]'); 
xlim([11 23.5]); ylim([-20.5 -9]);
set(gca,'xtick',10:1:55,'ytick',-38:1:-8,'tickdir','in','LineWidth',1.1,'FontSize',10,'FontName','Avenir'); 
grid on;

% print figure
%figurename=sprintf('type_cases/C5/quiver_movie_2023_200s_250430/positions_quiver_sZEROm_2023_200s_250430_%d.png',k+10);
figurename=sprintf('../../paperfigs/suppfig6.4_posstack_2023_200s_260130_%d.png',k+10);
print(gcf,'-dpng','-r300',figurename);  

%% ffmpeg line
%% ffmpeg -r 6 -f image2 -pattern_type glob -i "positions*.png" -vcodec libx264 -crf 25 -vf "crop=trunc(iw/2)*2:trunc(ih/2)*2" -pix_fmt yuv420p C6-2023-quiver.mp4

end