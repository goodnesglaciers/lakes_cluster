%% c_dot vs. epsilon_dot_lon vs. runoff over sections of the array
% 2024-01-23 LAS: plotting c_dot, longitudinal strain rates, runoff; need
% to fix en-dashes (LOL); need to understand error levels
% 2024-01-24 LAS: added en-dashes, location map, and c_dot labels; need to
% understand error levels!!!!
% 2024-03-21 LAS: updated after found error in daily_strain_rates calculations
% LAS 2024-03-26: tight and loose constraints for sliding least-squares,
% making a choice on which (station,days) for using which contraints
% LAS 2024-08-20: noting that this is the preferred 5-panel one for paperfig,
% removed small map
% LAS 2025-03-05: lots of small edits for paper version
% LAS 2025-04-21: update FASTER lake locations for map, 30-min vels
% LAS 2025-06-04: ice-flow directions on map-view panel
% LAS 25-06-12: +annotations
% LAS 26-01-16: +error envelopes on longitudinal strain rates; repository

close all; clear all

%% load runoff data
load('../../nevis_lakes_cluster/racmo_station_2023_300m_index.mat') % runoff_2023_nevis(:,ID(6,1))./10 = cm w.e. at SQ13
load('../../nevis_lakes_cluster/runoff_2023_nevis_noSK_300m.mat');
runoff_2023_nevis = runoff_2023_nevis.runoff;
% load station names
load station_names_2023.mat
num_sta = length(station_names); % number of stations

% load c_dot timeseries -- loose
load('daily_epsilon_zz_2023R_30min_BF2_UP4_sZERO_OUT6_clean_w36_t12_260119.mat')
daily_epsilon_zz_loose = daily_epsilon_zz;

% load c_dot timeseries -- tight
load('daily_epsilon_zz_2023R_30min_BF2_UP4_sZERO_OUT6_clean_w18_t6_260119.mat')
daily_epsilon_zz_tight = daily_epsilon_zz;

% load strain rates -- loose
load('daily_strain_rates_2023R_30min_BF2_UP4_sZERO_OUT6_clean_w36_t12_260119.mat')
daily_strain_rates_loose = daily_strain_rates_2023;
strain_time = daily_strain_rates_2023.time;
station_names_short = upper(daily_strain_rates_2023.station_names); % less SQ33

% load strain rates -- tight
load('daily_strain_rates_2023R_30min_BF2_UP4_sZERO_OUT6_clean_w18_t6_260119.mat')
daily_strain_rates_tight = daily_strain_rates_2023; 

% Tight constraints at station-days 2023:
tight_100s = 188:196;
tight_MLOW = 189:205;
tight_200s = 195:205;
tight_300s = 195:205;

%% COMBINE TIGHT AND LOOSE RECORDS FOR CDOT
% make timeseries logical for where to use 'tight' strain rates

for i=1:3 % MHIH, MLOW, QIET
    % first fill with loose constraints (longer sliding LS window)
    % loose daily_epsilon_zz
    daily_epsilon_zz(i).t22 = daily_epsilon_zz_loose(i).t22; % time
    daily_epsilon_zz(i).c_dot_delta_t_combo = daily_epsilon_zz_loose(i).c_dot_delta_t; % c_dot

    for j=1:1:length(daily_epsilon_zz(i).t22)
        if daily_epsilon_zz_tight(i).t22(j,1) >= tight_MLOW(1) && daily_epsilon_zz_tight(i).t22(j,1) <= tight_MLOW(end)+1
        % assign outputs from tight constraint (shorter sliding LS window)
        daily_epsilon_zz(i).c_dot_delta_t_combo(j,1) = daily_epsilon_zz_tight(i).c_dot_delta_t(j,1); % c_dot tight
        logical_combo(i).TF(j,1) = 1; % yes, use tight constraint
        else % keep loose constraints
        logical_combo(i).TF(j,1) = 0; % no, use loose constraint 
        end
    end
end

for i=4:9 % 100s
    % first fill with loose constraints
    % loose daily_epsilon_zz
    daily_epsilon_zz(i).t22 = daily_epsilon_zz_loose(i).t22; % time
    daily_epsilon_zz(i).c_dot_delta_t_combo = daily_epsilon_zz_loose(i).c_dot_delta_t; % c_dot
    
    for j=1:1:length(daily_epsilon_zz(i).t22)
        if daily_epsilon_zz_tight(i).t22(j,1) >= tight_100s(1) && daily_epsilon_zz_tight(i).t22(j,1) <= tight_100s(end)+1
        % assign outputs from tight constraint
        daily_epsilon_zz(i).c_dot_delta_t_combo(j,1) = daily_epsilon_zz_tight(i).c_dot_delta_t(j,1); % c_dot tight
        logical_combo(i).TF(j,1) = 1; % yes, use tight constraint
        else % keep loose constraints
            logical_combo(i).TF(j,1) = 0; % no, use loose constraint 
        end
    end
end

for i=10:16 % 200s
    % first fill with loose constraints
    % loose daily_epsilon_zz
    daily_epsilon_zz(i).t22 = daily_epsilon_zz_loose(i).t22; % time
    daily_epsilon_zz(i).c_dot_delta_t_combo = daily_epsilon_zz_loose(i).c_dot_delta_t; % c_dot
    
    for j=1:1:length(daily_epsilon_zz(i).t22)
        if daily_epsilon_zz_tight(i).t22(j,1) >= tight_200s(1) && daily_epsilon_zz_tight(i).t22(j,1) <= tight_200s(end)+1
        % assign outputs from tight constraint
        daily_epsilon_zz(i).c_dot_delta_t_combo(j,1) = daily_epsilon_zz_tight(i).c_dot_delta_t(j,1); % c_dot tight
        logical_combo(i).TF(j,1) = 1; % yes, use tight constraint
        else % keep loose constraints
            logical_combo(i).TF(j,1) = 0; % no, use loose constraint 
        end
    end
end

for i=17:num_sta-1 % 300s
    % first fill with loose constraints
    % loose daily_epsilon_zz
    daily_epsilon_zz(i).t22 = daily_epsilon_zz_loose(i).t22; % time
    daily_epsilon_zz(i).c_dot_delta_t_combo = daily_epsilon_zz_loose(i).c_dot_delta_t; % c_dot
    
    for j=1:1:length(daily_epsilon_zz(i).t22)
        if daily_epsilon_zz_tight(i).t22(j,1) >= tight_300s(1) && daily_epsilon_zz_tight(i).t22(j,1) <= tight_300s(end)+1
        % assign outputs from tight constraint
        daily_epsilon_zz(i).c_dot_delta_t_combo(j,1) = daily_epsilon_zz_tight(i).c_dot_delta_t(j,1); % c_dot tight
        logical_combo(i).TF(j,1) = 1; % yes, use tight constraint
        else % keep loose constraints
            logical_combo(i).TF(j,1) = 0; % no, use loose constraint 
        end
    end
end

%% COMBINE TIGHT AND LOOSE RECORDS FOR LONGITUDINAL STRAINS
% first start with loose constraints
daily_strain_rates_combo = daily_strain_rates_loose; 
% then fill in tight constraints
for i=1:1:num_sta-1 % first station
    for j=1:1:num_sta-1 % station's pair  
        for k=1:1:length(daily_strain_rates_combo.time) % over time
         % assign tight constraint if time of first station denoted as tight constraint
         if logical_combo(i).TF(k,1) == 1 
             % assign imputs from tight constraint   
             daily_strain_rates_combo.lon_yr_100s(i,j,k) = daily_strain_rates_tight.lon_yr_100s(i,j,k); % lons
             daily_strain_rates_combo.lon_yr_200s(i,j,k) = daily_strain_rates_tight.lon_yr_200s(i,j,k);
             daily_strain_rates_combo.lon_yr_300s(i,j,k) = daily_strain_rates_tight.lon_yr_300s(i,j,k);

             daily_strain_rates_combo.delta_lon_yr_100s(i,j,k) = daily_strain_rates_tight.delta_lon_yr_100s(i,j,k); % lons errors
             daily_strain_rates_combo.delta_lon_yr_200s(i,j,k) = daily_strain_rates_tight.delta_lon_yr_200s(i,j,k);
             daily_strain_rates_combo.delta_lon_yr_300s(i,j,k) = daily_strain_rates_tight.delta_lon_yr_300s(i,j,k);

             daily_strain_rates_combo.trans_yr_100s(i,j,k) = daily_strain_rates_tight.trans_yr_100s(i,j,k); % trans
             daily_strain_rates_combo.trans_yr_200s(i,j,k) = daily_strain_rates_tight.trans_yr_200s(i,j,k);
             daily_strain_rates_combo.trans_yr_300s(i,j,k) = daily_strain_rates_tight.trans_yr_300s(i,j,k);

             daily_strain_rates_combo.shear_yr_100s(i,j,k) = daily_strain_rates_tight.shear_yr_100s(i,j,k); % shears
             daily_strain_rates_combo.shear_yr_200s(i,j,k) = daily_strain_rates_tight.shear_yr_200s(i,j,k);
             daily_strain_rates_combo.shear_yr_300s(i,j,k) = daily_strain_rates_tight.shear_yr_300s(i,j,k);
        else % keep loose constraints
        end
        end
    end
end

% summer-average longitudinal strain rate error
for i=1:1:22 % first station
    for j=1:1:22 % station's pair 
    xxx1(i,j) = squeeze(nanmean(daily_strain_rates_combo.delta_lon_yr_100s(i,j,:)));
    xxx2(i,j) = squeeze(nanmean(daily_strain_rates_combo.delta_lon_yr_200s(i,j,:)));
    xxx3(i,j) = squeeze(nanmean(daily_strain_rates_combo.delta_lon_yr_300s(i,j,:)));
    end
end

%% load north lake geographic files and morlighem bed relative to North Lake
% load stations in polarstereographic coords
load polarstereo_stations_2023_short.mat
% load bedmap
load('../strain_rates_2022/BMv5_for_nevis_catchment.mat'); % origin = M1 moulin 
load('../strain_rates_2022/cmapland2.mat'); % bed topo colormap  
%  2023 FASTER boundaries and correct dates:
load('environs_lakes_2023B_S1S2MO_251003.mat') %  boundaries and correct dates
environs_lakes_2023B = environs_lakes; 
% load lake boundaries, centre points, drainage dates from FASTER for 2022
load('../strain_rates_2022/environs_lakes_2022B_250416.mat') % drainage dates & boundaries
environs_lakes_2022 = environs_lakes; 
% ROI that's relevant to the plotted region:
ROI_relevant_xv = [-15 -15 67 67]; % km in nevis-model space
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
%figure(2); plot(ROI_relevant_xv,ROI_relevant_yv)
%hold on; plot(environs_lakes_2023B.X_km,environs_lakes_2023B.Y_km,'o'); 

% XY relative to moulin M1 in Stevens et al. (2015)
origin = [68.72, -49.53]; % M1 moulin
% polarstereo conversion needed values
radius=6378137.0;    eccen=0.08181919;    lat_true=70;    lon_posy=-45;
[moulin_x,moulin_y] = polarstereo_fwd(origin(1),origin(2),radius,eccen,lat_true,lon_posy); % [ m ]
moulin_x_km = moulin_x./1e3; moulin_y_km = moulin_y./1e3; % UTM in KM of moulin
% load Nat/Laura's 2023 moulins
NT_23 = shaperead('../strain_rates_2022/NatalieTurner_Moulins_2022/moulins_2023.shp');
[moulins23_x,moulins23_y] = polarstereo_fwd([NT_23.Y],[NT_23.X],radius,eccen,lat_true,lon_posy); % [ m ]
moulins_x_2023 = (moulins23_x-moulin_x)./1e3; moulins_y_2023 = (moulins23_y-moulin_y)./1e3;
% load Nat/Laura's connections to lakes and moulins
load('../strain_rates_2022/NatalieTurner_Moulins_2022/laketypeing_2023_251002_result.mat');

%% plotting preliminaries 
lake_dates = [191.2, 196.75, 202.75, 200.50, 195.25]; % 100s, MHIH (IDs 150&152), 200s, 300s, IDs 154&162
lake_fill_dates = [179, 190, 183, 190, 190]; % 100s, MHIH, 200s, 300s, IDs 154&162
% cluster dates
dates_cluster = [191, 194.5, 200.5]; % C4, C5.N/C5.S, C6
dates_cluster_width = [3, 7, 1]; % C4, C5.N/C5.S, C6
% strain-rate time
delta_t = daily_epsilon_zz(1).delta_t;
% racmo
racmo_time = 1.5:1:365.5; 
% t_0 for strain_time
[ID_t0] = find(strain_time>=165.01, 1, 'first');

%% figure -- 2023 stack: maps, strain rate, runoff
figure(1); clf; 
set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',1.19*[0 0 18.3 24.7]);
strains_plot_height = 0.19;
axe1 = axes('Position',[0.027 0.5975 0.965 0.4015],'Box','on','NextPlot','add','XTickLabels',[],'YTickLabels',[]); % map
axe2 = axes('Position',[0.067 0.035+0.015+(strains_plot_height*(5/6))+(strains_plot_height*(4/6))+(strains_plot_height*(2/6)) 0.88 strains_plot_height*(5/6)],...
    'Box','on','NextPlot','add','XTickLabels',[]); % 100s
axeT = axes('Position',[0.067 0.035+0.010+(strains_plot_height*(5/6))+(strains_plot_height*(4/6)) 0.88 strains_plot_height*(2/6)],...
    'Box','on','NextPlot','add','XTickLabels',[]); % ties
axe3 = axes('Position',[0.067 0.035+0.005+(strains_plot_height*(5/6)) 0.88 strains_plot_height*(4/6)],'Box','on','NextPlot','add','XTickLabels',[]); % 200s
axe4 = axes('Position',[0.067 0.035 0.88 strains_plot_height*(5/6)],'Box','on','NextPlot','add'); % 300s

t_plot = [180 215]; % 2023
%t_plot = [180 220]; % 2022
% 2023
y_plot_strain = [-0.10 0.20]; % strain, runoff only
y_plot_strain_tick = 0.05;
y_plot_runoff = [0 20];
y_plot_runoff_tick = 4.0; 

% figure preliminaries
fontsize = 9; 
TriangleSize = 7;
DiamondSize = 0.5;
linewidth_strain = 1.3;
% dock at eel pond colormap :)
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

%%%%%%%%%%%%%%%%%%%%%%%%% 100s (in-basin) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(axe2) 
colororder({'k','k'})
yyaxis right
colororder('default'); CMap = get(gca, 'ColorOrder'); hold on
% runoff
patch_vert = nanmean(runoff_2023_nevis(:,ID(4:9,1)),2)./10; % RACMO in cm w.e.
patch_x1 = racmo_time-0.5+0.25; patch_x2 = patch_x1+0.5;
patch_y1 = ones(length(racmo_time),1).*y_plot_runoff(1);
patch_y2 = patch_y1 + patch_vert; 
for i=1:1:length(racmo_time)
    patch([patch_x1(i) patch_x2(i) patch_x2(i) patch_x1(i)],...
          [patch_y1(i) patch_y1(i) patch_y2(i) patch_y2(i)],...
          [0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.30); hold on;
end
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_runoff(1) y_plot_runoff(2)])
ylabel('Runoff at 950s [cm w.e.]','FontName','Avenir','FontSize',fontsize)
set(gca,'ytick',y_plot_runoff(1)+y_plot_runoff_tick:y_plot_runoff_tick:y_plot_runoff(2));
text(t_plot(2)-0.5,18,'within-basin 950s','FontName','Avenir','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','right')

text(dates_cluster(1)+1.1, y_plot_runoff(2)-1.4, 'C4','FontName','AvenirBold','FontSize',fontsize+1,...
    'FontAngle','italic','HorizontalAlignment','left','Rotation',0); 
text(dates_cluster(2)+2.6, y_plot_runoff(2)-1.4, 'C5','FontName','AvenirBold','FontSize',fontsize+1,...
    'FontAngle','italic','HorizontalAlignment','left','Rotation',0); 
text(dates_cluster(3)+0.1, y_plot_runoff(2)-1.4, 'C6','FontName','AvenirBold','FontSize',fontsize+1,...
    'FontAngle','italic','HorizontalAlignment','left','Rotation',0);

yyaxis left % \dot{epsilon}_{lon}
for i=1 % lakes filling and draining
    rectangle('Position',[lake_fill_dates(i) y_plot_strain(1) lake_dates(i)-lake_fill_dates(i) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',lake_blue,'EdgeColor','none'); 
    rectangle('Position',[dates_cluster(i) y_plot_strain(1) dates_cluster_width(i) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[1 1 0.7],'EdgeColor','none');
     text(180.25,-0.025,'950s filling...',...
         'FontName','Avenir','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','left')
end
for i=2:1:length(dates_cluster) % other drainages
    rectangle('Position',[dates_cluster(i) y_plot_strain(1) dates_cluster_width(i) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');
end

% strain rates
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(4,5,:)-daily_strain_rates_combo.lon_yr_100s(4,5,ID_t0)),...
    '-','LineWidth',linewidth_strain); hold on; %SQ11-12
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(6,7,:)-daily_strain_rates_combo.lon_yr_100s(6,7,ID_t0)),...
    '-','LineWidth',linewidth_strain); %SQ13-14
plot(strain_time(1:end), squeeze(daily_strain_rates_combo.lon_yr_100s(8,9,1:end)-daily_strain_rates_combo.lon_yr_100s(8,9,ID_t0)),...
    '-','LineWidth',linewidth_strain); %SQ15-16
% strain-rate error envelopes (any data gaps --> error estimate == 0)
tt = [strain_time; flipud(strain_time)];
ee = [squeeze(daily_strain_rates_combo.lon_yr_100s(4,5,:)+3.*daily_strain_rates_combo.delta_lon_yr_100s(4,5,:)-daily_strain_rates_combo.lon_yr_100s(4,5,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_100s(4,5,:)-3.*daily_strain_rates_combo.delta_lon_yr_100s(4,5,:)-daily_strain_rates_combo.lon_yr_100s(4,5,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_100s(6,7,:)+3.*daily_strain_rates_combo.delta_lon_yr_100s(6,7,:)-daily_strain_rates_combo.lon_yr_100s(6,7,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_100s(6,7,:)-3.*daily_strain_rates_combo.delta_lon_yr_100s(6,7,:)-daily_strain_rates_combo.lon_yr_100s(6,7,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_100s(8,9,:)+3.*daily_strain_rates_combo.delta_lon_yr_100s(8,9,:)-daily_strain_rates_combo.lon_yr_100s(8,9,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_100s(8,9,:)-3.*daily_strain_rates_combo.delta_lon_yr_100s(8,9,:)-daily_strain_rates_combo.lon_yr_100s(8,9,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
% strain rates
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(4,5,:)-daily_strain_rates_combo.lon_yr_100s(4,5,ID_t0)),'-','LineWidth',linewidth_strain,'Color',CMap(1,:)); hold on; %SQ11-12
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(6,7,:)-daily_strain_rates_combo.lon_yr_100s(6,7,ID_t0)),'-','LineWidth',linewidth_strain,'Color',CMap(2,:)); %SQ12-14
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(8,9,:)-daily_strain_rates_combo.lon_yr_100s(8,9,ID_t0)),'-','LineWidth',linewidth_strain,'Color',CMap(3,:)); %SQ15-16

plot([dates_cluster(3) dates_cluster(3)], [y_plot_strain(1) y_plot_strain(2)],':k','LineWidth',0.9); 
text(192, 0.125, {'L1A, L1B, & L1C';'HF drainages'},'FontName','Avenir','FontSize',fontsize,'FontAngle','italic','HorizontalAlignment','left')

colororder('default');
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_strain(1)+0.05 y_plot_strain(2)])
ylabel('$\dot{\epsilon}_{lon}(t)$  [ yr$^{-1}$ ]','Interpreter','latex','FontSize',fontsize+2)
legend(['SQ11' char(8211) '12'],['SQ13' char(8211) '14'],['SQ15' char(8211) '16'],'NumColumns',1,...
    'Location','NorthWest','EdgeColor','none','Color','none')
    % 'Position',[0.205 0.805 0.01 0.01],'EdgeColor','none')
grid on; grid minor; 
set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.6,'FontSize',fontsize,'MinorGridLineStyle',':','YMinorGrid','off'); 
set(gca,'ytick',y_plot_strain(1)+y_plot_strain_tick+y_plot_strain_tick:y_plot_strain_tick:y_plot_strain(2))
axe2.Layer = 'top';
text(177.45, 0.19, 'b','FontWeight','bold','FontSize',fontsize+3);

%%%%%%%%%%%%%%%%%%%%%%%%% TIEPOINTS (out-basin) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
axes(axeT) 
colororder({'k','k'})
yyaxis right 
colororder('default'); hold on;
% runoff
patch_vert = nanmean(runoff_2023_nevis(:,ID(2:3,1)),2)./10; % RACMO cm w.e.
patch_x1 = racmo_time-0.5+0.25; patch_x2 = patch_x1+0.5;
patch_y1 = ones(length(racmo_time),1).*y_plot_runoff(1);
patch_y2 = patch_y1 + patch_vert; 
for i=1:1:length(racmo_time)
    patch([patch_x1(i) patch_x2(i) patch_x2(i) patch_x1(i)],...
          [patch_y1(i) patch_y1(i) patch_y2(i) patch_y2(i)],...
          [0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.30); hold on;
end
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_runoff(1) 8])
ylabel('1050s','FontSize',fontsize,'FontName','Avenir')
set(gca,'ytick',y_plot_runoff(1)+y_plot_runoff_tick:y_plot_runoff_tick:8);
text(t_plot(2)-0.5,6,'tiepoints','FontName','Avenir','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','right')
        
yyaxis left % \dot{epsilon}_{lon}
for i=4 % lake filling and draining
    rectangle('Position',[dates_cluster(i-2) y_plot_strain(1) dates_cluster_width(i-2) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); hold on
end
for i=1:2:length(dates_cluster) % other drainages
    rectangle('Position',[dates_cluster(i) y_plot_strain(1) dates_cluster_width(i) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');
end

% Ties
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(3,5,:)-daily_strain_rates_combo.lon_yr_100s(3,5,ID_t0)),'-','LineWidth',linewidth_strain); hold on; %SQ12-QIET
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(3,14,:)-daily_strain_rates_combo.lon_yr_100s(3,14,ID_t0)),'-','LineWidth',linewidth_strain); %QIET–SQ25
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(2,10,:)-daily_strain_rates_combo.lon_yr_100s(2,10,ID_t0)),'-','LineWidth',linewidth_strain); %MLOW–SQ21
% strain-rate error envelopes (any data gaps --> error estimate == 0)
tt = [strain_time; flipud(strain_time)];
ee = [squeeze(daily_strain_rates_combo.lon_yr_100s(3,5,:)+3.*daily_strain_rates_combo.delta_lon_yr_100s(3,5,:)-daily_strain_rates_combo.lon_yr_100s(3,5,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_100s(3,5,:)-3.*daily_strain_rates_combo.delta_lon_yr_100s(3,5,:)-daily_strain_rates_combo.lon_yr_100s(3,5,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_100s(3,14,:)+3.*daily_strain_rates_combo.delta_lon_yr_100s(3,14,:)-daily_strain_rates_combo.lon_yr_100s(3,14,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_100s(3,14,:)-3.*daily_strain_rates_combo.delta_lon_yr_100s(3,14,:)-daily_strain_rates_combo.lon_yr_100s(3,14,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_100s(2,10,:)+3.*daily_strain_rates_combo.delta_lon_yr_100s(2,10,:)-daily_strain_rates_combo.lon_yr_100s(2,10,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_100s(2,10,:)-3.*daily_strain_rates_combo.delta_lon_yr_100s(2,10,:)-daily_strain_rates_combo.lon_yr_100s(2,10,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
% Ties
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(3,5,:)-daily_strain_rates_combo.lon_yr_100s(3,5,ID_t0)),'-','LineWidth',linewidth_strain,'Color',CMap(1,:)); hold on; %SQ12-QIET
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(3,14,:)-daily_strain_rates_combo.lon_yr_100s(3,14,ID_t0)),'-','LineWidth',linewidth_strain,'Color',CMap(2,:)); %QIET–SQ25
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(2,10,:)-daily_strain_rates_combo.lon_yr_100s(2,10,ID_t0)),'-','LineWidth',linewidth_strain,'Color',CMap(3,:)); %MLOW–SQ21

plot([dates_cluster(3) dates_cluster(3)], [y_plot_strain(1) y_plot_strain(2)],':k','LineWidth',0.9); 

colororder('default');
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_strain(1)+0.05 0.05])
ylabel('$\dot{\epsilon}_{lon}(t)$  [ yr$^{-1}$ ]','Interpreter','latex','FontSize',fontsize+1)
legend(['SQ12' char(8211) 'QIET'],['QIET' char(8211) 'SQ25'],['MLOW' char(8211) 'SQ21'],...
    'Location','NorthWest','Position',[0.147 0.365 0.01 0.01],'NumColumns',1,'EdgeColor','none','Color','none')
grid on; grid minor;
set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.6,'FontSize',fontsize,'MinorGridLineStyle',':','YMinorGrid','off'); 
set(gca,'ytick',y_plot_strain(1)+y_plot_strain_tick+y_plot_strain_tick:y_plot_strain_tick:y_plot_strain(2))
grid on; box on; axeT.Layer = 'top';
text(177.45, 0.07, 'c','FontWeight','bold','FontSize',fontsize+3);

%%%%%%%%%%%%%%%%%%%%%%%%% 200s (in-basin) %%%%%%%%%%%%%%%%%%%%%%%%%
axes(axe3) 
colororder({'k','k'})
yyaxis right 
colororder('default'); hold on;
% runoff
patch_vert = nanmean(runoff_2023_nevis(:,ID(10:16,1)),2)./10; % runoff in cm w.e.
patch_x1 = racmo_time-0.5+0.25; patch_x2 = patch_x1+0.5;
patch_y1 = ones(length(racmo_time),1).*y_plot_runoff(1);
patch_y2 = patch_y1 + patch_vert; 
for i=1:1:length(racmo_time)
    patch([patch_x1(i) patch_x2(i) patch_x2(i) patch_x1(i)],...
          [patch_y1(i) patch_y1(i) patch_y2(i) patch_y2(i)],...
          [0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.30); hold on;
end
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_runoff(1) 16])
ylabel('1150s [cm w.e.]','FontSize',fontsize,'FontName','Avenir')
set(gca,'ytick',y_plot_runoff(1)+y_plot_runoff_tick:y_plot_runoff_tick:16);
text(t_plot(2)-0.5,14,'within-basin 1150s','FontName','Avenir','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','right')

yyaxis left % \dot{epsilon}_{lon}
for i=2 % lakes filling and draining
    rectangle('Position',[lake_fill_dates(i+1) y_plot_strain(1) lake_dates(i+1)-lake_fill_dates(i+1) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',lake_blue,'EdgeColor','none'); % 200s
%     rectangle('Position',[dates_cluster(i) y_plot_strain(1) dates_cluster_width(i) y_plot_strain(2)-y_plot_strain(1)],...
%         'FaceColor',[1 1 0.7],'EdgeColor','none'); hold on
end
for i=1:1:length(dates_cluster) % other drainages
    rectangle('Position',[dates_cluster(i) y_plot_strain(1) dates_cluster_width(i) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');
end
text(lake_fill_dates(3)+0.25,-0.025,'1150s filling...',...
     'FontName','Avenir','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','left')

% 200s
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_200s(10,12,:)-daily_strain_rates_combo.lon_yr_200s(10,12,ID_t0)),...
    '-','LineWidth',linewidth_strain); hold on; %SQ21-23
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_200s(11,12,:)-daily_strain_rates_combo.lon_yr_200s(11,12,ID_t0)),...
    '-','LineWidth',linewidth_strain); %SQ22-23
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_200s(14,15,:)-daily_strain_rates_combo.lon_yr_200s(14,15,ID_t0)),...
    '-','LineWidth',linewidth_strain); %SQ25-26
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_200s(14,16,:)-daily_strain_rates_combo.lon_yr_200s(14,16,ID_t0)),...
    '-','LineWidth',linewidth_strain); hold on; %SQ25-27
% strain errors
tt = [strain_time; flipud(strain_time)];
ee = [squeeze(daily_strain_rates_combo.lon_yr_200s(10,12,:)+3.*daily_strain_rates_combo.delta_lon_yr_200s(10,12,:)-daily_strain_rates_combo.lon_yr_200s(10,12,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_200s(10,12,:)-3.*daily_strain_rates_combo.delta_lon_yr_200s(10,12,:)-daily_strain_rates_combo.lon_yr_200s(10,12,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_200s(11,12,:)+3.*daily_strain_rates_combo.delta_lon_yr_200s(11,12,:)-daily_strain_rates_combo.lon_yr_200s(11,12,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_200s(11,12,:)-3.*daily_strain_rates_combo.delta_lon_yr_200s(11,12,:)-daily_strain_rates_combo.lon_yr_200s(11,12,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_200s(14,15,:)+3.*daily_strain_rates_combo.delta_lon_yr_200s(14,15,:)-daily_strain_rates_combo.lon_yr_200s(14,15,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_200s(14,15,:)-3.*daily_strain_rates_combo.delta_lon_yr_200s(14,15,:)-daily_strain_rates_combo.lon_yr_200s(14,15,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_200s(14,16,:)+3.*daily_strain_rates_combo.delta_lon_yr_200s(14,16,:)-daily_strain_rates_combo.lon_yr_200s(14,16,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_200s(14,16,:)-3.*daily_strain_rates_combo.delta_lon_yr_200s(14,16,:)-daily_strain_rates_combo.lon_yr_200s(14,16,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
% 200s
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_200s(10,12,:)-daily_strain_rates_combo.lon_yr_200s(10,12,ID_t0)),...
    '-','LineWidth',linewidth_strain,'Color',CMap(1,:)); hold on; %SQ21-23
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_200s(11,12,:)-daily_strain_rates_combo.lon_yr_200s(11,12,ID_t0)),...
    '-','LineWidth',linewidth_strain,'Color',CMap(2,:)); %SQ22-23
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_200s(14,15,:)-daily_strain_rates_combo.lon_yr_200s(14,15,ID_t0)),...
    '-','LineWidth',linewidth_strain,'Color',CMap(3,:)); %SQ25-26
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_200s(14,16,:)-daily_strain_rates_combo.lon_yr_200s(14,16,ID_t0)),...
    '-','LineWidth',linewidth_strain,'Color',CMap(4,:)); hold on; %SQ25-27

plot([dates_cluster(3) dates_cluster(3)], [y_plot_strain(1) y_plot_strain(2)],':k','LineWidth',0.9); 

text(204, 0.05, {'L2A & L2B';'HF drainages'},'FontName','Avenir','FontSize',fontsize,'FontAngle','italic','HorizontalAlignment','left')

colororder('default');
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_strain(1) 0.1])
ylabel('$\dot{\epsilon}_{lon}(t)$  [ yr$^{-1}$ ]','Interpreter','latex','FontSize',fontsize+2)
legend(['SQ21' char(8211) '23'],['SQ22' char(8211) '23'],['SQ25' char(8211) '26'],['SQ25' char(8211) '27'],'NumColumns',2,...
    'Location','NorthWest','EdgeColor','none','Color','none')
    % 'Position',[0.2525 0.445 0.01 0.01],'EdgeColor','none') 
grid on; grid minor; 
set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.6,'FontSize',fontsize,'MinorGridLineStyle',':','YMinorGrid','off'); 
set(gca,'ytick',y_plot_strain(1)+y_plot_strain_tick:y_plot_strain_tick:0.1)
axe3.Layer = 'top';
text(177.45, 0.085, 'd','FontWeight','bold','FontSize',fontsize+3);

%%%%%%%%%%%%%%%%%%%%%%%%% 300s (in-basin) %%%%%%%%%%%%%%%%%%%%%%%%% 
axes(axe4) 
colororder({'k','k'})
yyaxis right 
colororder('default'); hold on;
% runoff
patch_vert = nanmean(runoff_2023_nevis(:,ID(17:23,1)),2)./10; % RACMO cm w.e.
patch_x1 = racmo_time-0.5+0.25; patch_x2 = patch_x1+0.5;
patch_y1 = ones(length(racmo_time),1).*y_plot_runoff(1);
patch_y2 = patch_y1 + patch_vert; 
for i=1:1:length(racmo_time)
    patch([patch_x1(i) patch_x2(i) patch_x2(i) patch_x1(i)],...
          [patch_y1(i) patch_y1(i) patch_y2(i) patch_y2(i)],...
          [0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.30); hold on;
end
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_runoff(1) y_plot_runoff(2)])
set(gca,'ytick',y_plot_runoff(1):y_plot_runoff_tick:y_plot_runoff(2));
ylabel('Runoff at 1350s [cm w.e.]','FontSize',fontsize,'FontName','Avenir')
text(t_plot(2)-0.5,18,'within-basin 1350s','FontName','Avenir','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','right')
        
yyaxis left % \dot{epsilon}_{lon}
for i=4 % lake filling and draining
    rectangle('Position',[lake_fill_dates(i) y_plot_strain(1) lake_dates(i)-lake_fill_dates(i) y_plot_strain(2)+0.1],...
        'FaceColor',lake_blue,'EdgeColor','none');
    rectangle('Position',[dates_cluster(i-2) y_plot_strain(1) dates_cluster_width(i-2) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[1 1 0.7],'EdgeColor','none'); hold on
end
for i=1:2 % other drainages
    rectangle('Position',[dates_cluster(i) y_plot_strain(1) dates_cluster_width(i) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');
end
rectangle('Position',[dates_cluster(3) y_plot_strain(1) dates_cluster_width(3) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[1 1 0.7],'EdgeColor','none'); hold on
text(lake_fill_dates(4)+0.25,-0.025,'1350s filling...',...
         'FontName','Avenir','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','left')

% 300s
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_300s(17,22,:)-daily_strain_rates_combo.lon_yr_300s(17,22,ID_t0)),...
    '-','LineWidth',linewidth_strain,'Color',CMap(1,:)); hold on; %SQ31-37
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_300s(18,19,:)-daily_strain_rates_combo.lon_yr_300s(18,19,ID_t0)),...
    '-','LineWidth',linewidth_strain,'Color',CMap(2,:)); hold on; %SQ32-34
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_300s(19,22,:)-daily_strain_rates_combo.lon_yr_300s(19,22,ID_t0)),...
    '-','LineWidth',linewidth_strain,'Color',CMap(3,:)); %SQ34-37
[index_182] = find((strain_time >= 182),1);
plot(strain_time(index_182:end), squeeze(daily_strain_rates_combo.lon_yr_300s(20,21,index_182:end)-daily_strain_rates_combo.lon_yr_300s(20,21,ID_t0)),...
    '-','LineWidth',linewidth_strain,'Color',CMap(4,:)); %SQ35-36
% strain errors
tt = [strain_time; flipud(strain_time)];
ee = [squeeze(daily_strain_rates_combo.lon_yr_300s(17,22,:)+3.*daily_strain_rates_combo.delta_lon_yr_300s(17,22,:)-daily_strain_rates_combo.lon_yr_300s(17,22,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_300s(17,22,:)-3.*daily_strain_rates_combo.delta_lon_yr_300s(17,22,:)-daily_strain_rates_combo.lon_yr_300s(17,22,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.75); hold on;
tt = [strain_time(index_182:end); flipud(strain_time(index_182:end))];
ee = [squeeze(daily_strain_rates_combo.lon_yr_300s(18,19,index_182:end)+3.*daily_strain_rates_combo.delta_lon_yr_300s(18,19,index_182:end)-daily_strain_rates_combo.lon_yr_300s(18,19,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_300s(18,19,index_182:end)-3.*daily_strain_rates_combo.delta_lon_yr_300s(18,19,index_182:end)-daily_strain_rates_combo.lon_yr_300s(18,19,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.75); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_300s(19,22,index_182:end)+3.*daily_strain_rates_combo.delta_lon_yr_300s(19,22,index_182:end)-daily_strain_rates_combo.lon_yr_300s(19,22,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_300s(19,22,index_182:end)-3.*daily_strain_rates_combo.delta_lon_yr_300s(19,22,index_182:end)-daily_strain_rates_combo.lon_yr_300s(19,22,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.75); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_300s(20,21,index_182:end)+3.*daily_strain_rates_combo.delta_lon_yr_300s(20,21,index_182:end)-daily_strain_rates_combo.lon_yr_300s(20,21,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_300s(20,21,index_182:end)-3.*daily_strain_rates_combo.delta_lon_yr_300s(20,21,index_182:end)-daily_strain_rates_combo.lon_yr_300s(20,21,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.75); hold on;
% 300s
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_300s(17,22,:)-daily_strain_rates_combo.lon_yr_300s(17,22,ID_t0)),...
    '-','LineWidth',linewidth_strain,'Color',CMap(1,:)); hold on; %SQ31-37
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_300s(18,19,:)-daily_strain_rates_combo.lon_yr_300s(18,19,ID_t0)),...
    '-','LineWidth',linewidth_strain,'Color',CMap(2,:)); hold on; %SQ32-34
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_300s(19,22,:)-daily_strain_rates_combo.lon_yr_300s(19,22,ID_t0)),...
    '-','LineWidth',linewidth_strain,'Color',CMap(3,:)); %SQ34-37
plot(strain_time(index_182:end), squeeze(daily_strain_rates_combo.lon_yr_300s(20,21,index_182:end)-daily_strain_rates_combo.lon_yr_300s(20,21,ID_t0)),...
    '-','LineWidth',linewidth_strain,'Color',CMap(4,:)); %SQ35-36

text(201.75, 0.072, {'L3A & L3B';'HF drainages'},'FontName','Avenir','FontSize',fontsize,'FontAngle','italic','HorizontalAlignment','left')

colororder('default');
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_strain(1) y_plot_strain(2)-y_plot_strain_tick])
ylabel('$\dot{\epsilon}_{lon}(t)$  [ yr$^{-1}$ ]','Interpreter','latex','FontSize',fontsize+2)
legend(['SQ31' char(8211) '37'],['SQ32' char(8211) '34'],['SQ34' char(8211) '37'],['SQ35' char(8211) '36'],'Location','NorthWest','NumColumns',1,...
            'Location','NorthWest','EdgeColor','none','color','none')
grid on; grid minor; 
set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.6,'FontSize',fontsize,'MinorGridLineStyle',':','YMinorGrid','off'); 
set(gca,'ytick',y_plot_strain(1):y_plot_strain_tick:y_plot_strain(2)-y_plot_strain_tick)
xlabel('Day of Year, 2023  [ UTC ]'); 
box on; axe4.Layer = 'top';
text(177.45, 0.15, 'e','FontWeight','bold','FontSize',fontsize+3);

%%%%%%%%%%%%%%%%%%%%%%%%%% MAP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% full ROI
plot_x = [-15 67];
plot_y = [-44 7];

axes(axe1) % 2023
%xlabel('East-West [ km ]','FontSize',fontsize);
%ylabel('North-South [ km ]','FontSize',fontsize);
set(gca,'xaxislocation','top','FontSize',fontsize)

% ice-sheet bed contour map, colorbar
surf(BMv5_for_nevis_catchment.X_km/1e3,BMv5_for_nevis_catchment.Y_km./1e3,...
   BMv5_for_nevis_catchment.B_km,'EdgeColor','None','Facealpha',0.5); % plot surface below station points
view([0 90])
caxis([-600 200]); colormap(bone)

% t2=colorbar('east');
% hold all
% set(get(t2,'xlabel'),'String',{'Bed Elevation  [ m a.s.l. ]'},'FontSize',fontsize,'Fontname','Avenir');
% set(t2, 'Position', [0.97 0.85 .006 0.14],'YTick',-600:200:200,...
%    'YTickLabel',[-600,-400,-200,0,200],'TickDirection','both');

hold all
% surface ice-sheet contours
surf_contours = 0:100:1900; % m a.s.l.
[C, h] = contour(BMv5_for_nevis_catchment.X_km./1e3, BMv5_for_nevis_catchment.Y_km./1e3, ...
    BMv5_for_nevis_catchment.S_km, surf_contours,'Color',[0.7 0.7 0.7],'LineWidth',0.8);
clabel(C,h,'FontSize',7,'Color',[0.7 0.7 0.7],'FontName','Avenir','LabelSpacing',400)

% % ice-flow direction
% span = 6; tick_scale = 1.5;
% quiver3(stresses_nlake_2022Annual_S12.XB_strain_grid_sub(1:span:end,1:span:end)./1e3,...
%     stresses_nlake_2022Annual_S12.YB_strain_grid_sub(1:span:end,1:span:end)./1e3,...
%     1e3.*ones(size(stresses_nlake_2022Annual_S12.YB_strain_grid_sub(1:span:end,1:span:end))),...
%     tick_scale.*cosd(stresses_nlake_2022Annual_S12.flow_angle(1:span:end,1:span:end)),...
%     tick_scale.*sind(stresses_nlake_2022Annual_S12.flow_angle(1:span:end,1:span:end)),...
%     0.*sind(stresses_nlake_2022Annual_S12.flow_angle(1:span:end,1:span:end)),...
%     'LineWidth',1.1,'ShowArrowHead','off','Color',[0.6 0.2 0.0],'AutoScale','off');

% array GNSS
plot3(xy_sta_23_short(:,1),xy_sta_23_short(:,2),3000.*ones(num_sta-1,1),'k^',...
    'MarkerSize',TriangleSize-2,'MarkerFaceColor','w','markerEdgeColor',metal,...
    'LineWidth',1.2); 


% ROI that's relevant to the GNSS array = [-30 100] in X; [-40 10] in Y
% if the lake location is not within relevant catchment, throw it out !! 
ROI_relevant_xv = [-15 60 104 104 15 -10 -33 -33 -15]; % km in nevis-model space
ROI_relevant_yv = [15 4 -8 -50 -50 -41 -33 -1 15];
plot3(ROI_relevant_xv,ROI_relevant_yv,1e4.*[1 1 1 1 1 1 1 1 1],'-','LineWidth',2,'Color',[0.6 0.2 0]);

% lakes: FASTER S2 2023
[HF_lakes_to_plot] = find(environs_lakes_2023B.laketypeing_dates(:,7)==1); % just HF lakes  
[MU_lakes_to_plot] = find(environs_lakes_2023B.laketypeing_dates(:,7)==2); % just moulin-draining lakes 
[lakes_to_plot] = find(~isnan(environs_lakes_2023B.laketypeing_dates(:,7))); % all good lakes   
[NE_lakes_to_plot] = find(environs_lakes_2023B.laketypeing_dates(:,7)==5); % just no-exit freezers

% plotting in order for legend
fill3(environs_lakes_2023B.boundaries(lakes_to_plot(1)).XY_km_local(:,1), ...
        environs_lakes_2023B.boundaries(lakes_to_plot(1)).XY_km_local(:,2), ...
        3000.*ones(length(environs_lakes_2023B.boundaries(lakes_to_plot(1)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',lake_blue(1:3));
fill3(environs_lakes_2023B.boundaries(HF_lakes_to_plot(1)).XY_km_local(:,1), ...
        environs_lakes_2023B.boundaries(HF_lakes_to_plot(1)).XY_km_local(:,2), ...
        3000.*ones(length(environs_lakes_2023B.boundaries(HF_lakes_to_plot(1)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',goldenrod);
fill3(environs_lakes_2023B.boundaries(MU_lakes_to_plot(1)).XY_km_local(:,1), ...
        environs_lakes_2023B.boundaries(MU_lakes_to_plot(1)).XY_km_local(:,2), ...
        3000.*ones(length(environs_lakes_2023B.boundaries(MU_lakes_to_plot(1)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',[0.4 0.8 0.8]);
fill3(environs_lakes_2023B.boundaries(NE_lakes_to_plot(1)).XY_km_local(:,1), ...
        environs_lakes_2023B.boundaries(NE_lakes_to_plot(1)).XY_km_local(:,2), ...
        3000.*ones(length(environs_lakes_2023B.boundaries(NE_lakes_to_plot(1)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',plum(1:3));
% Natalie Turner's lake connections for overspill lakes, JUST IN ROI
for i=43
    if strcmp(result_23{i,3}, 'L')
    plot3([environs_lakes_2022.X_km(i) environs_lakes_2022.X_km(result_23{i,4})],...
          [environs_lakes_2022.Y_km(i) environs_lakes_2022.Y_km(result_23{i,4})],[300 300],...
          '-','LineWidth',0.6,'Color',[0 0.55 1])
    else end
end

% cluster cOneXiOnS
grey_connections = lake_blue(1:3).*1.3; % BLUE (MAIN TEXT)
% 10s
cluster_10s = [44, 66, 75, 81, 86, 103, 131, 134];
cc1 = plot3([environs_lakes_2023B.X_km(cluster_10s(1)) environs_lakes_2023B.X_km(cluster_10s(3))],...
    [environs_lakes_2023B.Y_km(cluster_10s(1)) environs_lakes_2023B.Y_km(cluster_10s(3))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
cc12 =plot3([environs_lakes_2023B.X_km(cluster_10s(2)) 0],...
    [environs_lakes_2023B.Y_km(cluster_10s(2)) 4],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
cc13 =plot3([environs_lakes_2023B.X_km(cluster_10s(3)) environs_lakes_2023B.X_km(cluster_10s(4))],...
    [environs_lakes_2023B.Y_km(cluster_10s(3)) environs_lakes_2023B.Y_km(cluster_10s(4))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
cc13 =plot3([environs_lakes_2023B.X_km(cluster_10s(5)) environs_lakes_2023B.X_km(cluster_10s(4))],...
    [environs_lakes_2023B.Y_km(cluster_10s(5)) environs_lakes_2023B.Y_km(cluster_10s(4))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
cc13 =plot3([environs_lakes_2023B.X_km(cluster_10s(3)) environs_lakes_2023B.X_km(cluster_10s(6))],...
    [environs_lakes_2023B.Y_km(cluster_10s(3)) environs_lakes_2023B.Y_km(cluster_10s(6))],[2e3 2e3],...
    ':','LineWidth',2,'Color',grey_connections);
cc13 =plot3([environs_lakes_2023B.X_km(cluster_10s(5)) environs_lakes_2023B.X_km(cluster_10s(7))],...
    [environs_lakes_2023B.Y_km(cluster_10s(5)) environs_lakes_2023B.Y_km(cluster_10s(7))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
cc13 =plot3([environs_lakes_2023B.X_km(cluster_10s(5)) 0.1],...
    [environs_lakes_2023B.Y_km(cluster_10s(5)) 4.0],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections); %L1C
cc13 =plot3([environs_lakes_2023B.X_km(cluster_10s(8)) environs_lakes_2023B.X_km(cluster_10s(7))],...
    [environs_lakes_2023B.Y_km(cluster_10s(8)) environs_lakes_2023B.Y_km(cluster_10s(7))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
text(-6, -5, 5e3, 'C4','FontSize',fontsize+7,'FontName','AvenirBold','FontAngle','italic',...
    'FontWeight','bold','Color',grey_connections);
% #mid
cluster_mid = [146, 150, 152, 154, 162, 171, 179, 181, 194, 196];
cc1 = plot3([environs_lakes_2023B.X_km(cluster_mid(1)) environs_lakes_2023B.X_km(cluster_mid(6))],...
    [environs_lakes_2023B.Y_km(cluster_mid(1)) environs_lakes_2023B.Y_km(cluster_mid(6))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
cc1 = plot3([environs_lakes_2023B.X_km(cluster_mid(8)) environs_lakes_2023B.X_km(cluster_mid(6))],...
    [environs_lakes_2023B.Y_km(cluster_mid(8)) environs_lakes_2023B.Y_km(cluster_mid(6))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
cc11 =plot3([environs_lakes_2023B.X_km(cluster_mid(2)) environs_lakes_2023B.X_km(cluster_mid(3))],...
    [environs_lakes_2023B.Y_km(cluster_mid(2)) environs_lakes_2023B.Y_km(cluster_mid(3))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
cc12 =plot3([environs_lakes_2023B.X_km(cluster_mid(3)) environs_lakes_2023B.X_km(cluster_mid(8))],...
    [environs_lakes_2023B.Y_km(cluster_mid(3)) environs_lakes_2023B.Y_km(cluster_mid(8))],[2e3 2e3],...
    ':','LineWidth',2,'Color',grey_connections);
cc13 =plot3([environs_lakes_2023B.X_km(cluster_mid(4)) environs_lakes_2023B.X_km(cluster_mid(7))],...
    [environs_lakes_2023B.Y_km(cluster_mid(4)) environs_lakes_2023B.Y_km(cluster_mid(7))],[2e3 2e3],...
    ':','LineWidth',2,'Color',grey_connections); %% TAKEN OUT OF C5.N when checking dates!
cc13 =plot3([environs_lakes_2023B.X_km(cluster_mid(9)) environs_lakes_2023B.X_km(cluster_mid(5))],...
    [environs_lakes_2023B.Y_km(cluster_mid(9)) environs_lakes_2023B.Y_km(cluster_mid(5))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
cc13 =plot3([environs_lakes_2023B.X_km(cluster_mid(9)) environs_lakes_2023B.X_km(cluster_mid(10))],...
    [environs_lakes_2023B.Y_km(cluster_mid(9)) environs_lakes_2023B.Y_km(cluster_mid(10))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
cc13 =plot3([environs_lakes_2023B.X_km(cluster_mid(7)) environs_lakes_2023B.X_km(cluster_mid(10))],...
    [environs_lakes_2023B.Y_km(cluster_mid(7)) environs_lakes_2023B.Y_km(cluster_mid(10))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
text(23, -37, 5e3,'C5.S','FontSize',fontsize+7,'FontName','AvenirBold','FontAngle','italic',...
    'FontWeight','bold','Color',grey_connections);
text(30, -1.85, 5e3, 'C5.N','FontSize',fontsize+7,'FontName','AvenirBold','FontAngle','italic',...
    'FontWeight','bold','Color',grey_connections);
% 30s
cluster_30s = [208, 201, 228, 238, 273, 275];
cc2 = plot3([environs_lakes_2023B.X_km(cluster_30s(1)) environs_lakes_2023B.X_km(cluster_30s(2))],...
    [environs_lakes_2023B.Y_km(cluster_30s(1)) environs_lakes_2023B.Y_km(cluster_30s(2))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
cc21 =plot3([environs_lakes_2023B.X_km(cluster_30s(2)) environs_lakes_2023B.X_km(cluster_30s(3))],...
    [environs_lakes_2023B.Y_km(cluster_30s(2)) environs_lakes_2023B.Y_km(cluster_30s(3))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
cc22 =plot3([environs_lakes_2023B.X_km(cluster_30s(3)) environs_lakes_2023B.X_km(cluster_30s(4))],...
    [environs_lakes_2023B.Y_km(cluster_30s(3)) environs_lakes_2023B.Y_km(cluster_30s(4))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
cc23 =plot3([environs_lakes_2023B.X_km(cluster_30s(1)) environs_lakes_2023B.X_km(cluster_30s(5))],...
    [environs_lakes_2023B.Y_km(cluster_30s(1)) environs_lakes_2023B.Y_km(cluster_30s(5))],[2e3 2e3],...
    ':','LineWidth',2,'Color',grey_connections);
cc24 =plot3([environs_lakes_2023B.X_km(cluster_30s(6)) environs_lakes_2023B.X_km(cluster_30s(5))],...
    [environs_lakes_2023B.Y_km(cluster_30s(6)) environs_lakes_2023B.Y_km(cluster_30s(5))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
text(45.5, -27,5e3, 'C6','FontSize',fontsize+7,'FontName','AvenirBold','FontAngle','italic',...
    'FontWeight','bold','Color',grey_connections);

% Natalie Turner's lake connections for overspill lakes, JUST IN ROI
for i=1:1:406
    if strcmp(result_23{i,3}, 'L')
    plot3([environs_lakes_2023B.X_km(i) environs_lakes_2023B.X_km(result_23{i,4})],...
          [environs_lakes_2023B.Y_km(i) environs_lakes_2023B.Y_km(result_23{i,4})],[300 300],...
          '-','LineWidth',0.4,'Color',[0 0.55 1])
    elseif strcmp(result_23{i,3}, 'M')
    plot3([environs_lakes_2023B.X_km(i) moulins_x_2023(result_23{i,4})],...
        [environs_lakes_2023B.Y_km(i) moulins_y_2023(result_23{i,4})],[300 300],...
        '-','LineWidth',0.4,'Color',[0 0.55 1])
    else end
end
% Natalie Turner's moulins, JUST IN ROI
plot3(moulins_x_2023, moulins_y_2023,300.*ones(length(moulins_x_2023),1),'ko',...
     'MarkerSize',TriangleSize-4,'MarkerFaceColor','w','markerEdgeColor','k',...
     'LineWidth',0.5);


% fill those lakes
for i=1:1:length(lakes_to_plot)
    fill3(environs_lakes_2023B.boundaries(lakes_to_plot(i)).XY_km_local(:,1), ...
        environs_lakes_2023B.boundaries(lakes_to_plot(i)).XY_km_local(:,2), ...
        3000.*ones(length(environs_lakes_2023B.boundaries(lakes_to_plot(i)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',lake_blue(1:3));
end
for i=1:1:length(HF_lakes_to_plot)
    fill3(environs_lakes_2023B.boundaries(HF_lakes_to_plot(i)).XY_km_local(:,1), ...
        environs_lakes_2023B.boundaries(HF_lakes_to_plot(i)).XY_km_local(:,2), ...
        3000.*ones(length(environs_lakes_2023B.boundaries(HF_lakes_to_plot(i)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',goldenrod);
end
for i=1:1:length(MU_lakes_to_plot)
    fill3(environs_lakes_2023B.boundaries(MU_lakes_to_plot(i)).XY_km_local(:,1), ...
        environs_lakes_2023B.boundaries(MU_lakes_to_plot(i)).XY_km_local(:,2), ...
        3000.*ones(length(environs_lakes_2023B.boundaries(MU_lakes_to_plot(i)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',[0.4 0.8 0.8]);
end
for i=1:1:length(NE_lakes_to_plot) % No-exit freezers 
    fill3(environs_lakes_2023B.boundaries(NE_lakes_to_plot(i)).XY_km_local(:,1), ...
        environs_lakes_2023B.boundaries(NE_lakes_to_plot(i)).XY_km_local(:,2), ...
        3000.*ones(length(environs_lakes_2023B.boundaries(NE_lakes_to_plot(i)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',plum(1:3));
end

% L1C from 2022 outline (2022 ID = 77)
fill3(environs_lakes_2022.boundaries(77).XY_km_local(:,1), ...
        environs_lakes_2022.boundaries(77).XY_km_local(:,2), ...
        3000.*ones(length(environs_lakes_2022.boundaries(77).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',goldenrod);
text(environs_lakes_2022.X_km(77)-1.0, environs_lakes_2022.Y_km(77),5e3,...
        sprintf('191'), ...
        'FontSize',fontsize-2,'FontName','Helvetica','FontWeight','bold')

% L3A from 2022 outline (2022 ID = 232)
fill3(environs_lakes_2022.boundaries(232).XY_km_local(:,1), ...
        environs_lakes_2022.boundaries(232).XY_km_local(:,2), ...
        3000.*ones(length(environs_lakes_2022.boundaries(232).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',goldenrod);

% HF dates of drainages range (plot as text)
for i=1:1:length(HF_lakes_to_plot)
    j=HF_lakes_to_plot(i);  
       if environs_lakes_2023B.laketypeing_dates(j,3)+1 == environs_lakes_2023B.laketypeing_dates(j,4)
      t=text(environs_lakes_2023B.X_km(j)-1.0, environs_lakes_2023B.Y_km(j),5e3,...
        sprintf('%d',environs_lakes_2023B.laketypeing_dates(j,3)+1), ...
        'FontSize',fontsize-2,'FontName','Helvetica','FontWeight','bold');
        t.Clipping = "on";
       else
      t=text(environs_lakes_2023B.X_km(j)-1.5, environs_lakes_2023B.Y_km(j),5e3,...
        sprintf('%d\x2013%d',environs_lakes_2023B.laketypeing_dates(j,3)+1,environs_lakes_2023B.laketypeing_dates(j,4)), ...
        'FontSize',fontsize-2,'FontName','Helvetica','FontWeight','bold');
      t.Clipping = "on";
   end
end

% Moulin dates of drainages (plot as text)
for i=1:1:length(MU_lakes_to_plot)
    j=MU_lakes_to_plot(i);
       if environs_lakes_2023B.laketypeing_dates(j,3)+1 == environs_lakes_2023B.laketypeing_dates(j,4)
      t=text(environs_lakes_2023B.X_km(j)-1.0, environs_lakes_2023B.Y_km(j),5e3,...
        sprintf('%d',environs_lakes_2023B.laketypeing_dates(j,3)+1), ...
        'FontSize',fontsize-2,'FontName','Helvetica','FontWeight','bold');
      t.Clipping = "on";
   else 
    t=text(environs_lakes_2023B.X_km(j)-1.5, environs_lakes_2023B.Y_km(j),5e3,...
        sprintf('%d\x2013%d',environs_lakes_2023B.laketypeing_dates(j,3)+1,environs_lakes_2023B.laketypeing_dates(j,4)), ...
        'FontSize',fontsize-2,'FontName','Helvetica','FontWeight','bold');
    t.Clipping = "on";
   end
end

% site names
for i=1; text(xy_sta_23_short(i,1)-3.75,xy_sta_23_short(i,2)+0.55,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=2; text(xy_sta_23_short(i,1)-4.25,xy_sta_23_short(i,2)+0.55,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=3; text(xy_sta_23_short(i,1)+1.25,xy_sta_23_short(i,2)+0.55,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=4; text(xy_sta_23_short(i,1)-3.0,xy_sta_23_short(i,2)-1.55,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic','Rotation',30); end
for i=5; text(xy_sta_23_short(i,1)+0.5,xy_sta_23_short(i,2)-1.00,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=6; text(xy_sta_23_short(i,1)+1.25,xy_sta_23_short(i,2)+0.10,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=9; text(xy_sta_23_short(i,1)-3.25,xy_sta_23_short(i,2)+1.00,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=7; text(xy_sta_23_short(i,1)-3.75,xy_sta_23_short(i,2)-0.1,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=8; text(xy_sta_23_short(i,1)-3.75,xy_sta_23_short(i,2)+0.75,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=10:11; text(xy_sta_23_short(i,1)-3.75,xy_sta_23_short(i,2)-0.1,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=12:13; text(xy_sta_23_short(i,1)+1.25,xy_sta_23_short(i,2)-0.55,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=14; text(xy_sta_23_short(i,1)-3.75,xy_sta_23_short(i,2)+0.55,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=15:16; text(xy_sta_23_short(i,1)+1.25,xy_sta_23_short(i,2)+0.55,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=17; text(xy_sta_23_short(i,1)-0.75,xy_sta_23_short(i,2)-1.25,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=18; text(xy_sta_23_short(i,1)-3.5,xy_sta_23_short(i,2)+0.5,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=19; text(xy_sta_23_short(i,1)+0.5,xy_sta_23_short(i,2)-1.00,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic','rotation',-15); end
for i=20; text(xy_sta_23_short(i,1)-4.00,xy_sta_23_short(i,2)-0.5,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=21; text(xy_sta_23_short(i,1)+0.75,xy_sta_23_short(i,2)+0.75,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=22; text(xy_sta_23_short(i,1)+1.00,xy_sta_23_short(i,2)+0.55,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic','Rotation',45); end
% lake names
text(2, 4 ,5e3, 'L1C','Fontname','AvenirBold','FontSize',fontsize+2,'FontAngle','italic','Color','k')
text(-2.5, 0.75 ,5e3, 'L1A','Fontname','AvenirBold','FontSize',fontsize+2,'FontAngle','italic','Color','k')
text(0.25, -4.5, 5e3, 'L1B','Fontname','AvenirBold','FontSize',fontsize+2,'FontAngle','italic','Color','k')
text(23, -11.5, 5e3, 'L2A','Fontname','AvenirBold','FontSize',fontsize+2,'FontAngle','italic','Color','k')
text(14.5, -20 ,5e3, 'L2B','Fontname','AvenirBold','FontSize',fontsize+2,'FontAngle','italic','Color','k')
text(25, -17.5 ,5e3, 'L2C','Fontname','AvenirBold','FontSize',fontsize+2,'FontAngle','italic','Color','k')
text(38, -21.25, 5e3, 'L3A','Fontname','AvenirBold','FontSize',fontsize+2,'FontAngle','italic','Color','k')
text(31.5, -27.5, 5e3, 'L3B','Fontname','AvenirBold','FontSize',fontsize+2,'FontAngle','italic','Color','k')

% scalebar
plot3([50 60],[1 1],[2e3 2e3],'k-','LineWidth',2);
text(53, 2, 5e3, '10 km','FontName','AvenirBold','FontSize',fontsize+1)

lgd = legend('Bed Elevation','Ice Elevation','GNSS','ROI',...
    'Overspill drainage','HF drainage','Moulin drainage','No-exit, frozen',...
    "Overspill connections",'HF Clusters ({\itp}-value<0.05)');
set(lgd,'Position',[0.51 0.5725 0.01 0.01],'FontSize',fontsize,'FontName','Avenir',...
    'NumColumns',5,'Box','off','TextColor',lake_blue(1:3))

%axis equal; 
grid on; box on;
text(-17, 6.5, 'a','FontWeight','bold','FontSize',fontsize+3);
text(62.4, 5.5, 10e3, '2023','FontSize',fontsize+4,'FontName','Helvetica','FontWeight','bold')
set(gca,'FontName','Avenir','FontSize',fontsize,'xtick',-80:10:110,...
    'ytick',-50:10:30,'LineWidth',0.6)
xlim([plot_x(1) plot_x(2)]); ylim([plot_y(1) plot_y(2)]);
axe1.Layer = 'top';   axe1.ClippingStyle = "rectangle";  

%% print figure
print(gcf,'-dpng','-r300',sprintf('../../paperfigs/paperfig6_distillations_2023_260116.png')); 
