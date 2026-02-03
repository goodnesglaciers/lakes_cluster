%% epsilon_dot_lon vs. runoff over sections of the array
% 2024-01-23 LAS: plotting c_dot, longitudinal strain rates, runoff; need
% to fix en-dashes (LOL); need to understand error levels
% 2024-01-24 LAS: added en-dashes, location map, and c_dot labels; need to
% understand error levels!!!!
% 2024-03-21 LAS: updated after found error in daily_strain_rates calculations
% LAS 2024-03-26: tight and loose constraints for sliding least-squares,
% making a choice on which (station,days) for using which contraints
% LAS 2024-08-20: noting that this is the preferred 5-panel one for paperfig,
% removed small map
% LAS 2024-11-22: map + 3 panels of strain rates
% LAS 25-04-21: lots of small betterments for paper, 30-min clean vels
% LAS 25-06-12: +annotations
% LAS 26-01-16: +error envelopes on longitudinal strain rates

close all; clear all

% load runoff data
load('../../nevis_lakes_cluster/racmo_station_2022_index.mat') % runoff_2022_nevis(:,ID(6,1))./10 = cm w.e. at SQ13
load('../../nevis_lakes_cluster/runoff_2022_nevis');
% load station names
load station_names.mat; station_names = upper(station_names);

% load c_dot timeseries -- loose
load('daily_epsilon_zz_2022R_BF2_30min_BF2_UP4_sZERO_clean_w36_t12_260119.mat') % 30-min
daily_epsilon_zz_loose = daily_epsilon_zz;

% load c_dot timeseries -- tight
load('daily_epsilon_zz_2022R_BF2_30min_BF2_UP4_sZERO_clean_w18_t6_260119.mat') % 30-min
daily_epsilon_zz_tight = daily_epsilon_zz;

% load strain rates -- loose
load('daily_strain_rates_2022R_30min_BF2_UP4_sZERO_clean_w36_t12_260119.mat') % 30-min
daily_strain_rates_loose = daily_strain_rates_2023;
strain_time = daily_strain_rates_loose.time;
station_names_short = upper(daily_strain_rates_loose.station_names); % less SQ33

% load strain rates -- tight
load('daily_strain_rates_2022R_30min_BF2_UP4_sZERO_clean_w18_t6_260119.mat') % 30-min
daily_strain_rates_tight = daily_strain_rates_2023;

% Tight constraints at specific station-days: 
tight_100s = 187:196;
tight_200s = 206:217;
tight_300s = 206:212;
tight_MLOW = 207:217; 
tight_MHIH = 207:217;
tight_QIET = 192:196;
tight_QIET2 = 210:218;

% Tight constraints on all days: (spurious diurnal signals due to
% multipath are not related to glacier physics/hydrology) 
% tight_100s = 170:230;
% tight_200s = 170:230;
% tight_300s = 170:230;
% tight_MLOW = 170:230; 
% tight_MHIH = 170:230;
% tight_QIET = 170:230;
% tight_QIET2 = 170:230;

%% COMBINE TIGHT AND LOOSE RECORDS FOR CDOT
% make timeseries logical for where to use 'tight' strain rates
for i=1 % MHIH
    % first fill with loose constraints (longer sliding LS window)
    % loose daily_epsilon_zz
    daily_epsilon_zz(i).t22 = daily_epsilon_zz_loose(i).t22; % time
    daily_epsilon_zz(i).c_dot_delta_t_combo = daily_epsilon_zz_loose(i).c_dot_delta_t; % c_dot
    
    for j=1:1:length(daily_epsilon_zz(i).t22)
        if daily_epsilon_zz_tight(i).t22(j,1) >= tight_MHIH(1) && daily_epsilon_zz_tight(i).t22(j,1) <= tight_MHIH(end)+1
        % assign outputs from tight constraint (shorter sliding LS window)
        daily_epsilon_zz(i).c_dot_delta_t_combo(j,1) = daily_epsilon_zz_tight(i).c_dot_delta_t(j,1); % c_dot tight
        logical_combo(i).TF(j,1) = 1; % yes, use tight constraint

        else % keep loose constraints
        logical_combo(i).TF(j,1) = 0; % no, use loose constraint    
        end
    end
end

for i=2 % MLOW
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

for i=3 % QIET
    % first fill with loose constraints (longer sliding LS window)
    % loose daily_epsilon_zz
    daily_epsilon_zz(i).t22 = daily_epsilon_zz_loose(i).t22; % time
    daily_epsilon_zz(i).c_dot_delta_t_combo = daily_epsilon_zz_loose(i).c_dot_delta_t; % c_dot
    
    for j=1:1:length(daily_epsilon_zz(i).t22)
        if daily_epsilon_zz_tight(i).t22(j,1) >= tight_QIET(1) && daily_epsilon_zz_tight(i).t22(j,1) <= tight_QIET(end)+1
        % assign outputs from tight constraint (shorter sliding LS window)
        daily_epsilon_zz(i).c_dot_delta_t_combo(j,1) = daily_epsilon_zz_tight(i).c_dot_delta_t(j,1); % c_dot tight
        logical_combo(i).TF(j,1) = 1; % yes, use tight constraint
        elseif daily_epsilon_zz_tight(i).t22(j,1) >= tight_QIET2(1) && daily_epsilon_zz_tight(i).t22(j,1) <= tight_QIET2(end)+1
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

for i=17:21 % 300s
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
for i=1:1:21 % first station
    for j=1:1:21 % station's pair  
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

% summer-average longitudinal strain rate error (for table S10)
for i=1:1:21 % first station
    for j=1:1:21 % station's pair 
    xxx1(i,j) = squeeze(nanmean(daily_strain_rates_combo.delta_lon_yr_100s(i,j,:)));
    xxx2(i,j) = squeeze(nanmean(daily_strain_rates_combo.delta_lon_yr_200s(i,j,:)));
    xxx3(i,j) = squeeze(nanmean(daily_strain_rates_combo.delta_lon_yr_300s(i,j,:)));
    end
end

%% load north lake geographic files and morlighem bed relative to North Lake
% origin = [68.72, -49.53];
load polarstereo_stations_2022_short.mat
% load bedmap
load BMv5_for_nevis_catchment.mat % origin = M1 moulin 
load cmapland2.mat % bed topo colormap
% load lake boundaries, centre points, drainage dates from FASTER for 2022
load('environs_lakes_2022B_250416.mat') % drainage dates & boundaries
environs_lakes_2022 = environs_lakes; 

% ROI that's relevant to the plotted region:
ROI_relevant_xv = [-27.75 -27.75 60 60]; % km in nevis-model space
ROI_relevant_yv = [-36.8 9 9 -36.8];
for i=1:1:length(environs_lakes_2022.H)
    if inpolygon(environs_lakes_2022.X_km(i),environs_lakes_2022.Y_km(i),...
            ROI_relevant_xv,ROI_relevant_yv) == 0
        % throw a NaN
        environs_lakes_2022.H(i) = NaN;
        environs_lakes_2022.S(i) = NaN;
        environs_lakes_2022.X_km(i) = NaN;
        environs_lakes_2022.Y_km(i) = NaN;
        environs_lakes_2022.lat(i) = NaN;
        environs_lakes_2022.lon(i) = NaN;
        environs_lakes_2022.drainage_type_num(i) = NaN;
        environs_lakes_2022.laketypeing_dates(i,1:7) = [NaN, NaN, NaN, NaN, NaN, NaN, NaN];
    else
    end
end
% redefine lake_typeing to just include lakes within nevis domain
lake_typeing = environs_lakes_2022.laketypeing_dates; % after throwing out lakes on the nevis domain border
% polygon ROI check
% figure(2); plot(ROI_relevant_xv,ROI_relevant_yv)
% hold on; plot(environs_lakes_2022.X_km,environs_lakes_2022.Y_km,'o'); 

% XY relative to moulin M1 in Stevens et al. (2015)
origin = [68.72, -49.53]; % M1 moulin
% polarstereo conversion needed values
radius=6378137.0;    eccen=0.08181919;    lat_true=70;    lon_posy=-45;
[moulin_x,moulin_y] = polarstereo_fwd(origin(1),origin(2),radius,eccen,lat_true,lon_posy); % [ m ]
moulin_x_km = moulin_x./1e3; moulin_y_km = moulin_y./1e3; % UTM in KM of moulin
% load Natalie Turner's moulins
NT = shaperead('NatalieTurner_Moulins_2022/moulins.shp');
[moulins_x,moulins_y] = polarstereo_fwd([NT.Y],[NT.X],radius,eccen,lat_true,lon_posy); % [ m ]
moulins_x_2022 = (moulins_x-moulin_x)./1e3; moulins_y_2022 = (moulins_y-moulin_y)./1e3;
% load Natalie Turner's connections to lakes and moulins
load('laketypeing_2022_250930_result.mat')

%% plotting preliminaries 
t_plot = [150 254];

lake_dates = [195.0, 210, 214.5, 210]; % 100s, MHIH, 200s, L3B
lake_fill_dates = [176, 184, 184, 190]; % 100s, MHIH, 200s, 300s

dates_cluster = [195, 207, 214]; % C1, C2.N/C2.S, C3
dates_cluster_width = [1, 4, 2]; % C1, C2.N/C2.S, C3

% SQ13 not past 200 (heightened diurnals as antenna tipped over)
[index_200] = find((daily_epsilon_zz(6).t22 >= 200),1);
daily_epsilon_zz(6).c_dot_delta_t_combo(index_200:end) = NaN;
daily_epsilon_zz(6).u_s(index_200:end) = NaN;
daily_epsilon_zz(6).w_s(index_200:end) = NaN;
daily_epsilon_zz(6).epsilon_dot_lon(index_200:end) = NaN;
daily_epsilon_zz(6).epsilon_dot_trans(index_200:end) = NaN;
daily_epsilon_zz(6).epsilon_dot_zz(index_200:end) = NaN;
% SQ15 not past 230 (fell over -- DIVA DOWN!)
[index_230] = find((daily_epsilon_zz(8).t22 >= 230),1);
daily_epsilon_zz(8).c_dot_delta_t_combo(index_230:end) = NaN;
daily_epsilon_zz(8).u_s(index_230:end) = NaN;
daily_epsilon_zz(8).w_s(index_230:end) = NaN;
daily_epsilon_zz(8).epsilon_dot_lon(index_230:end) = NaN;
daily_epsilon_zz(8).epsilon_dot_trans(index_230:end) = NaN;
daily_epsilon_zz(7).epsilon_dot_trans(index_230:end) = NaN; % SQ14 uses SQ15 for yy strain rates
daily_epsilon_zz(8).epsilon_dot_zz(index_230:end) = NaN;
daily_epsilon_zz(7).epsilon_dot_zz(index_230:end) = NaN; % SQ14 uses SQ15 for yy strain rates
% strain-rate time
delta_t = daily_epsilon_zz(1).delta_t;
% racmo
racmo_time = 1.5:1:334.5; 
% t_0 for strain_time (DOY 165)
[ID_t0] = find(strain_time>=165.01, 1, 'first');

%% figure -- 2022 stack: maps, strain rate, runoff
figure(1); clf; 
set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',1.19*[0 0 18.3 24.7]);
strains_plot_height = 0.1625;
axe1 = axes('Position',[0.027 0.6075 0.965 0.3915],'Box','on','NextPlot','add','XTickLabels',[],'YTickLabels',[]); % map
axe2 = axes('Position',[0.067 0.035+0.015+(strains_plot_height*(3/6))+(strains_plot_height*(3/6))+(strains_plot_height*(7/6)) 0.88 strains_plot_height],...
    'Box','on','NextPlot','add','XTickLabels',[]); % 100s
axeT = axes('Position',[0.067 0.035+0.010+(strains_plot_height*(3/6))+(strains_plot_height*(7/6)) 0.88 strains_plot_height*(3/6)],...
    'Box','on','NextPlot','add','XTickLabels',[]); % ties
axe3 = axes('Position',[0.067 0.035+0.005+(strains_plot_height*(3/6)) 0.88 strains_plot_height*(7/6)],'Box','on','NextPlot','add','XTickLabels',[]); % 200s
axe4 = axes('Position',[0.067 0.035 0.88 strains_plot_height*(3/6)],'Box','on','NextPlot','add'); % 300s
%
t_plot = [180 225]; % meltseason 
% t_plot = [180 250]; % meltseason + late-summer melt event
y_plot_strain = [-0.06 0.12]; % strain, runoff only
y_plot_strain_tick = 0.03;
y_plot_runoff = [0 18];
y_plot_runoff_tick = 3.0; 

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
patch_vert = nanmean(runoff_2022_nevis(:,ID(4:9,1)),2)./10; % RACMO in cm w.e.
patch_x1 = racmo_time-0.5+0.25; patch_x2 = patch_x1+0.5;
patch_y1 = ones(length(racmo_time),1).*y_plot_runoff(1);
patch_y2 = patch_y1 + patch_vert; 
for i=1:1:length(racmo_time)
    patch([patch_x1(i) patch_x2(i) patch_x2(i) patch_x1(i)],...
          [patch_y1(i) patch_y1(i) patch_y2(i) patch_y2(i)],...
          [0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.30); hold on;
end
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_runoff(1) y_plot_runoff(2)])
ylabel('950s [cm w.e.]','FontName','Avenir','FontSize',fontsize)
set(gca,'ytick',y_plot_runoff(1)+y_plot_runoff_tick:y_plot_runoff_tick:y_plot_runoff(2));
text(t_plot(2)-0.5,16.5,'within-basin 950s','FontName','Avenir','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','right')

text(dates_cluster(1)+1.2, y_plot_runoff(2)-1.4, 'C1','FontName','AvenirBold','FontSize',fontsize+1,...
    'FontAngle','italic','HorizontalAlignment','left','Rotation',0); 
text(dates_cluster(2)+1.5, y_plot_runoff(2)-1.4, 'C2','FontName','AvenirBold','FontSize',fontsize+1,...
    'FontAngle','italic','HorizontalAlignment','left','Rotation',0); 
text(dates_cluster(3)+0.5, y_plot_runoff(2)-1.4, 'C3','FontName','AvenirBold','FontSize',fontsize+1,...
    'FontAngle','italic','HorizontalAlignment','left','Rotation',0);

yyaxis left % \dot{epsilon}_{lon}
for i=1 % lakes filling and draining
    rectangle('Position',[lake_fill_dates(i) y_plot_strain(1) lake_dates(i)-lake_fill_dates(i) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',lake_blue,'EdgeColor','none'); 
    rectangle('Position',[dates_cluster(i) y_plot_strain(1) dates_cluster_width(i) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[1 1 0.7],'EdgeColor','none');
     text(180.4,-0.015,'950s filling...',...
         'FontName','Avenir','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','left')
end
for i=2:1:length(dates_cluster) % other drainages
    rectangle('Position',[dates_cluster(i) y_plot_strain(1) dates_cluster_width(i) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');
end

% strain rates
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(4,5,:)-daily_strain_rates_combo.lon_yr_100s(4,5,ID_t0)),'-','LineWidth',linewidth_strain,'Color',CMap(1,:)); hold on; %SQ11-12
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(5,7,:)-daily_strain_rates_combo.lon_yr_100s(5,7,ID_t0)),'-','LineWidth',linewidth_strain,'Color',CMap(2,:)); %SQ12-14
plot(strain_time(1:1950), squeeze(daily_strain_rates_combo.lon_yr_100s(8,9,1:1950)-daily_strain_rates_combo.lon_yr_100s(8,9,ID_t0)),'-','LineWidth',linewidth_strain,'Color',CMap(3,:)); %SQ15-16
% strain-rate error envelopes (any data gaps --> error estimate == 0)
tt = [strain_time; flipud(strain_time)];
ee = [squeeze(daily_strain_rates_combo.lon_yr_100s(4,5,:)+3.*daily_strain_rates_combo.delta_lon_yr_100s(4,5,:)-daily_strain_rates_combo.lon_yr_100s(4,5,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_100s(4,5,:)-3.*daily_strain_rates_combo.delta_lon_yr_100s(4,5,:)-daily_strain_rates_combo.lon_yr_100s(4,5,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_100s(5,7,:)+3.*daily_strain_rates_combo.delta_lon_yr_100s(5,7,:)-daily_strain_rates_combo.lon_yr_100s(5,7,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_100s(5,7,:)-3.*daily_strain_rates_combo.delta_lon_yr_100s(5,7,:)-daily_strain_rates_combo.lon_yr_100s(5,7,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
patch([strain_time(1:1950); flipud(strain_time(1:1950))], ...
    [squeeze(daily_strain_rates_combo.lon_yr_100s(8,9,1:1950)+3.*daily_strain_rates_combo.delta_lon_yr_100s(8,9,1:1950)-daily_strain_rates_combo.lon_yr_100s(8,9,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_100s(8,9,1:1950)-3.*daily_strain_rates_combo.delta_lon_yr_100s(8,9,1:1950)-daily_strain_rates_combo.lon_yr_100s(8,9,ID_t0)))],...
    [0.7 0.7 0.7],'EdgeColor', 'none','FaceAlpha', 0.75); 
% strain rates
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(4,5,:)-daily_strain_rates_combo.lon_yr_100s(4,5,ID_t0)),'-','LineWidth',linewidth_strain,'Color',CMap(1,:)); hold on; %SQ11-12
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(5,7,:)-daily_strain_rates_combo.lon_yr_100s(5,7,ID_t0)),'-','LineWidth',linewidth_strain,'Color',CMap(2,:)); %SQ12-14
plot(strain_time(1:2000), squeeze(daily_strain_rates_combo.lon_yr_100s(8,9,1:2000)-daily_strain_rates_combo.lon_yr_100s(8,9,ID_t0)),'-','LineWidth',linewidth_strain,'Color',CMap(3,:)); %SQ15-16

text(196.25, 0.045, {'L1A & L1B';'HF drainages'},'FontName','Avenir','FontSize',fontsize,'FontAngle','italic','HorizontalAlignment','left')

colororder('default');
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_strain(1) y_plot_strain(2)])
ylabel('$\dot{\epsilon}_{lon}(t)$  [ yr$^{-1}$ ]','Interpreter','latex','FontSize',fontsize+1)
legend(['SQ11' char(8211) '12'],['SQ12' char(8211) '14'],['SQ15' char(8211) '16'],'NumColumns',1,...
    'Location','NorthWest','EdgeColor','none','Color','none')
grid on; grid minor;
set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.6,'FontSize',fontsize,'MinorGridLineStyle',':','YMinorGrid','off'); 
set(gca,'ytick',y_plot_strain(1)+y_plot_strain_tick:y_plot_strain_tick:y_plot_strain(2))
grid on; box on; axe2.Layer = 'top';
text(176.75, 0.11, 'b','FontWeight','bold','FontSize',fontsize+3);

%%%%%%%%%%%%%%%%%%%%%%%%% TIEPOINTS (out-basin) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
axes(axeT) 
colororder({'k','k'})
yyaxis right 
colororder('default'); hold on;
% runoff
patch_vert = nanmean(runoff_2022_nevis(:,ID(2:3,1)),2)./10; % RACMO cm w.e.
patch_x1 = racmo_time-0.5+0.25; patch_x2 = patch_x1+0.5;
patch_y1 = ones(length(racmo_time),1).*y_plot_runoff(1);
patch_y2 = patch_y1 + patch_vert; 
for i=1:1:length(racmo_time)
    patch([patch_x1(i) patch_x2(i) patch_x2(i) patch_x1(i)],...
          [patch_y1(i) patch_y1(i) patch_y2(i) patch_y2(i)],...
          [0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.30); hold on;
end
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_runoff(1) 9])
ylabel('1050s [cm w.e.]','FontSize',fontsize,'FontName','Avenir')
set(gca,'ytick',y_plot_runoff(1)+y_plot_runoff_tick:y_plot_runoff_tick:9);
text(t_plot(2)-0.5,7.5,'tiepoints','FontName','Avenir','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','right')
        
yyaxis left % \dot{epsilon}_{lon}
for i=4 % lake filling and draining
    rectangle('Position',[dates_cluster(i-2) y_plot_strain(1) dates_cluster_width(i-2) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); hold on
end
for i=1:2:length(dates_cluster) % other drainages
    rectangle('Position',[dates_cluster(i) y_plot_strain(1) dates_cluster_width(i) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');
end

% Ties's not before 2022/160
[index_184] = find((strain_time >= 160),1);
plot(strain_time(index_184:end), squeeze(daily_strain_rates_combo.lon_yr_100s(3,5,index_184:end)-daily_strain_rates_combo.lon_yr_100s(3,5,ID_t0)),'-','LineWidth',linewidth_strain); hold on; %SQ12-QIET
plot(strain_time(index_184:end), squeeze(daily_strain_rates_combo.lon_yr_100s(3,14,index_184:end)-daily_strain_rates_combo.lon_yr_100s(3,14,ID_t0)),'-','LineWidth',linewidth_strain); %QIET–SQ25
plot(strain_time(index_184:end), squeeze(daily_strain_rates_combo.lon_yr_100s(2,10,index_184:end)-daily_strain_rates_combo.lon_yr_100s(2,10,ID_t0)),'-','LineWidth',linewidth_strain); %MLOW–SQ21
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
% Ties's not before 2022/160
plot(strain_time(index_184:end), squeeze(daily_strain_rates_combo.lon_yr_100s(3,5,index_184:end)-daily_strain_rates_combo.lon_yr_100s(3,5,ID_t0)),'-','LineWidth',linewidth_strain,'Color',CMap(1,:)); hold on; %SQ12-QIET
plot(strain_time(index_184:end), squeeze(daily_strain_rates_combo.lon_yr_100s(3,14,index_184:end)-daily_strain_rates_combo.lon_yr_100s(3,14,ID_t0)),'-','LineWidth',linewidth_strain,'Color',CMap(2,:)); %QIET–SQ25
plot(strain_time(index_184:end), squeeze(daily_strain_rates_combo.lon_yr_100s(2,10,index_184:end)-daily_strain_rates_combo.lon_yr_100s(2,10,ID_t0)),'-','LineWidth',linewidth_strain,'Color',CMap(3,:)); %MLOW–SQ21

colororder('default');
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_strain(1) 0.03])
ylabel('$\dot{\epsilon}_{lon}(t)$  [ yr$^{-1}$ ]','Interpreter','latex','FontSize',fontsize+1)
legend(['SQ12' char(8211) 'QIET'],['QIET' char(8211) 'SQ25'],['MLOW' char(8211) 'SQ21'],...
    'Location','SouthWest','NumColumns',1,'EdgeColor','none','Color','none')
grid on; grid minor;
set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.6,'FontSize',fontsize,'MinorGridLineStyle',':','YMinorGrid','off'); 
set(gca,'ytick',y_plot_strain(1)+y_plot_strain_tick:y_plot_strain_tick:y_plot_strain(2))
grid on; box on; axeT.Layer = 'top';
text(176.75, 0.04, 'c','FontWeight','bold','FontSize',fontsize+3);


%%%%%%%%%%%%%%%%%%%%%%%%% 200s (in-basin) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y_plot_strain = [-0.06 0.15]; % strain, runoff only
y_plot_strain_tick = 0.03;
y_plot_runoff = [0 21];
y_plot_runoff_tick = 3.0; 

axes(axe3) 
colororder({'k','k'})
yyaxis right 
colororder('default'); hold on;
% runoff
patch_vert = nanmean(runoff_2022_nevis(:,ID(10:16,1)),2)./10; % runoff in cm w.e.
patch_x1 = racmo_time-0.5+0.25; patch_x2 = patch_x1+0.5;
patch_y1 = ones(length(racmo_time),1).*y_plot_runoff(1);
patch_y2 = patch_y1 + patch_vert; 
for i=1:1:length(racmo_time)
    patch([patch_x1(i) patch_x2(i) patch_x2(i) patch_x1(i)],...
          [patch_y1(i) patch_y1(i) patch_y2(i) patch_y2(i)],...
          [0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.30); hold on;
end
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_runoff(1) 18])
ylabel('1150s [cm w.e.]','FontSize',fontsize,'FontName','Avenir')
set(gca,'ytick',y_plot_runoff(1)+y_plot_runoff_tick:y_plot_runoff_tick:18);
text(t_plot(2)-0.5,16.5,'within-basin 1150s','FontName','Avenir','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','right')

yyaxis left % \dot{epsilon}_{lon}
for i=3 % lakes filling and draining
    rectangle('Position',[lake_fill_dates(i) y_plot_strain(1) lake_dates(i)-lake_fill_dates(i) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',lake_blue,'EdgeColor','none');
    rectangle('Position',[dates_cluster(i) y_plot_strain(1) dates_cluster_width(i) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[1 1 0.7],'EdgeColor','none'); hold on
     text(lake_fill_dates(i)+0.5,+0.015,'1150s filling...',...
         'FontName','Avenir','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','left')
end
for i=1:2 % other drainages
    rectangle('Position',[dates_cluster(i) y_plot_strain(1) dates_cluster_width(i) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');
end

% 200's not before 2022/182
[index_182] = find((strain_time >= 182),1);
plot(strain_time(index_182:end), squeeze(daily_strain_rates_combo.lon_yr_200s(10,12,index_182:end)-daily_strain_rates_combo.lon_yr_200s(10,12,ID_t0)),'-','LineWidth',linewidth_strain); hold on; %SQ21-23
plot(strain_time(index_182:end), squeeze(daily_strain_rates_combo.lon_yr_200s(11,12,index_182:end)-daily_strain_rates_combo.lon_yr_200s(11,12,ID_t0)),'-','LineWidth',linewidth_strain); %SQ22-23
plot(strain_time(index_182:end), squeeze(daily_strain_rates_combo.lon_yr_200s(14,15,index_182:end)-daily_strain_rates_combo.lon_yr_200s(14,15,ID_t0)),'-','LineWidth',linewidth_strain); %SQ25-26
plot(strain_time(index_182:end), squeeze(daily_strain_rates_combo.lon_yr_200s(14,16,index_182:end)-daily_strain_rates_combo.lon_yr_200s(14,16,ID_t0)),'-','LineWidth',linewidth_strain); hold on; %SQ25-27
% strain-rate error envelopes (any data gaps --> error estimate == 0)
tt = [strain_time(index_182:end); flipud(strain_time(index_182:end))];
ee = [squeeze(daily_strain_rates_combo.lon_yr_200s(10,12,index_182:end)+3.*daily_strain_rates_combo.delta_lon_yr_200s(10,12,index_182:end)-daily_strain_rates_combo.lon_yr_200s(10,12,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_200s(10,12,index_182:end)-3.*daily_strain_rates_combo.delta_lon_yr_200s(10,12,index_182:end)-daily_strain_rates_combo.lon_yr_200s(10,12,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_200s(11,12,index_182:end)+3.*daily_strain_rates_combo.delta_lon_yr_200s(11,12,index_182:end)-daily_strain_rates_combo.lon_yr_200s(11,12,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_200s(11,12,index_182:end)-3.*daily_strain_rates_combo.delta_lon_yr_200s(11,12,index_182:end)-daily_strain_rates_combo.lon_yr_200s(11,12,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_200s(14,15,index_182:end)+3.*daily_strain_rates_combo.delta_lon_yr_200s(14,15,index_182:end)-daily_strain_rates_combo.lon_yr_200s(14,15,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_200s(14,15,index_182:end)-3.*daily_strain_rates_combo.delta_lon_yr_200s(14,15,index_182:end)-daily_strain_rates_combo.lon_yr_200s(14,15,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_200s(14,16,index_182:end)+3.*daily_strain_rates_combo.delta_lon_yr_200s(14,16,index_182:end)-daily_strain_rates_combo.lon_yr_200s(14,16,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_200s(14,16,index_182:end)-3.*daily_strain_rates_combo.delta_lon_yr_200s(14,16,index_182:end)-daily_strain_rates_combo.lon_yr_200s(14,16,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
% 200's not before 2022/182
plot(strain_time(index_182:end), squeeze(daily_strain_rates_combo.lon_yr_200s(10,12,index_182:end)-daily_strain_rates_combo.lon_yr_200s(10,12,ID_t0)),'-','LineWidth',linewidth_strain,'Color',CMap(1,:)); hold on; %SQ21-23
plot(strain_time(index_182:end), squeeze(daily_strain_rates_combo.lon_yr_200s(11,12,index_182:end)-daily_strain_rates_combo.lon_yr_200s(11,12,ID_t0)),'-','LineWidth',linewidth_strain,'Color',CMap(2,:)); %SQ22-23
plot(strain_time(index_182:end), squeeze(daily_strain_rates_combo.lon_yr_200s(14,15,index_182:end)-daily_strain_rates_combo.lon_yr_200s(14,15,ID_t0)),'-','LineWidth',linewidth_strain,'Color',CMap(3,:)); %SQ25-26
plot(strain_time(index_182:end), squeeze(daily_strain_rates_combo.lon_yr_200s(14,16,index_182:end)-daily_strain_rates_combo.lon_yr_200s(14,16,ID_t0)),'-','LineWidth',linewidth_strain,'Color',CMap(4,:)); hold on; %SQ25-27

text(216.25, 0.045, {'L2A & L2B';'HF drainages'},'FontName','Avenir','FontSize',fontsize,'FontAngle','italic','HorizontalAlignment','left')

colororder('default');
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_strain(1)+y_plot_strain_tick y_plot_strain(2)])
ylabel('$\dot{\epsilon}_{lon}(t)$  [ yr$^{-1}$ ]','Interpreter','latex','FontSize',fontsize+1)
legend(['SQ21' char(8211) '23'],['SQ22' char(8211) '23'],['SQ25' char(8211) '26'],['SQ25' char(8211) '27'],'NumColumns',1,...
    'Location','NorthWest','EdgeColor','none','Color','none')
grid on; grid minor;
set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.6,'FontSize',fontsize,'MinorGridLineStyle',':','YMinorGrid','off'); 
set(gca,'ytick',y_plot_strain(1)+y_plot_strain_tick+y_plot_strain_tick:y_plot_strain_tick:y_plot_strain(2))
grid on; box on; axe3.Layer = 'top';
text(176.75, 0.135, 'd','FontWeight','bold','FontSize',fontsize+3);

%%%%%%%%%%%%%%%%%%%%%%%%% 300s (in-basin) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
axes(axe4) 
colororder({'k','k'})
yyaxis right 
colororder('default'); hold on;
% runoff
patch_vert = nanmean(runoff_2022_nevis(:,ID(17:21,1)),2)./10; % RACMO cm w.e.
patch_x1 = racmo_time-0.5+0.25; patch_x2 = patch_x1+0.5;
patch_y1 = ones(length(racmo_time),1).*y_plot_runoff(1);
patch_y2 = patch_y1 + patch_vert; 
for i=1:1:length(racmo_time)
    patch([patch_x1(i) patch_x2(i) patch_x2(i) patch_x1(i)],...
          [patch_y1(i) patch_y1(i) patch_y2(i) patch_y2(i)],...
          [0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.30); hold on;
end
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_runoff(1) 9])
ylabel('Runoff at 1350s [cm w.e.]','FontSize',fontsize,'FontName','Avenir')
set(gca,'ytick',y_plot_runoff(1):y_plot_runoff_tick:9);
text(t_plot(2)-0.5,7.5,'within-basin 1350s','FontName','Avenir','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','right')
        
yyaxis left % \dot{epsilon}_{lon}
for i=4 % lake filling and draining
    rectangle('Position',[lake_fill_dates(i) y_plot_strain(1) lake_dates(i)-lake_fill_dates(i) 0.06-y_plot_strain(1)],...
        'FaceColor',lake_blue,'EdgeColor','none');
    rectangle('Position',[dates_cluster(i-2) y_plot_strain(1) dates_cluster_width(i-2) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[1 1 0.7],'EdgeColor','none'); hold on
end
for i=1:2:length(dates_cluster) % other drainages
    rectangle('Position',[dates_cluster(i) y_plot_strain(1) dates_cluster_width(i) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[0.95 0.95 0.95],'EdgeColor','none');
end
     text(lake_fill_dates(4)+0.5,+0.015,'1350s filling...',...
         'FontName','Avenir','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','left')

% 300's not before 2022/184
[index_184] = find((strain_time >= 184),1);
plot(strain_time(index_184:end), squeeze(daily_strain_rates_combo.lon_yr_300s(17,18,index_184:end)-daily_strain_rates_combo.lon_yr_300s(17,18,ID_t0)),'-','LineWidth',linewidth_strain,'Color',CMap(1,:)); hold on; %SQ31-34
plot(strain_time(index_184:end), squeeze(daily_strain_rates_combo.lon_yr_300s(19,20,index_184:end)-daily_strain_rates_combo.lon_yr_300s(19,20,ID_t0)),'-','LineWidth',linewidth_strain,'Color',CMap(2,:)); %SQ35-36
plot(strain_time(index_184:end), squeeze(daily_strain_rates_combo.lon_yr_300s(19,21,index_184:end)-daily_strain_rates_combo.lon_yr_300s(19,21,ID_t0)),'-','LineWidth',linewidth_strain,'Color',CMap(3,:)); %SQ35-37
% strain-rate error envelopes (any data gaps --> error estimate == 0)
tt = [strain_time(index_184:end); flipud(strain_time(index_184:end))];
ee = [squeeze(daily_strain_rates_combo.lon_yr_300s(17,18,index_184:end)+3.*daily_strain_rates_combo.delta_lon_yr_300s(17,18,index_184:end)-daily_strain_rates_combo.lon_yr_300s(17,18,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_300s(17,18,index_184:end)-3.*daily_strain_rates_combo.delta_lon_yr_300s(17,18,index_184:end)-daily_strain_rates_combo.lon_yr_300s(17,18,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.75); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_300s(19,20,index_184:end)+3.*daily_strain_rates_combo.delta_lon_yr_300s(19,20,index_184:end)-daily_strain_rates_combo.lon_yr_300s(19,20,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_300s(19,20,index_184:end)-3.*daily_strain_rates_combo.delta_lon_yr_300s(19,20,index_184:end)-daily_strain_rates_combo.lon_yr_300s(19,20,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.75); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_300s(19,21,index_184:end)+3.*daily_strain_rates_combo.delta_lon_yr_300s(19,21,index_184:end)-daily_strain_rates_combo.lon_yr_300s(19,21,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_300s(19,21,index_184:end)-3.*daily_strain_rates_combo.delta_lon_yr_300s(19,21,index_184:end)-daily_strain_rates_combo.lon_yr_300s(19,21,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.75); hold on;
% 300's not before 2022/184
plot(strain_time(index_184:end), squeeze(daily_strain_rates_combo.lon_yr_300s(17,18,index_184:end)-daily_strain_rates_combo.lon_yr_300s(17,18,ID_t0)),'-','LineWidth',linewidth_strain,'Color',CMap(1,:)); hold on; %SQ31-34
plot(strain_time(index_184:end), squeeze(daily_strain_rates_combo.lon_yr_300s(19,20,index_184:end)-daily_strain_rates_combo.lon_yr_300s(19,20,ID_t0)),'-','LineWidth',linewidth_strain,'Color',CMap(2,:)); %SQ35-36
plot(strain_time(index_184:end), squeeze(daily_strain_rates_combo.lon_yr_300s(19,21,index_184:end)-daily_strain_rates_combo.lon_yr_300s(19,21,ID_t0)),'-','LineWidth',linewidth_strain,'Color',CMap(3,:)); %SQ35-37

text(211.25, 0.029, {'L3B moulin';'drainage'},'FontName','Avenir','FontSize',fontsize,'FontAngle','italic','HorizontalAlignment','left')

colororder('default');
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_strain(1)+0.03 0.06])
ylabel('$\dot{\epsilon}_{lon}(t)$  [ yr$^{-1}$ ]','Interpreter','latex','FontSize',fontsize+1)
legend(['SQ31' char(8211) '34'],'SQ35-36','SQ35-37','Location','NorthWest','NumColumns',1,...
    'Location','NorthWest','EdgeColor','none','Color','none')
grid on; grid minor;
set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.6,'FontSize',fontsize,'MinorGridLineStyle',':','YMinorGrid','off'); 
set(gca,'ytick',-0.03:y_plot_strain_tick:y_plot_strain(2))
xlabel('Day of Year, 2022  [ UTC ]'); 
box on; axe4.Layer = 'top';
text(176.75, 0.07, 'e','FontWeight','bold','FontSize',fontsize+3);

%%%%%%%%%%%%%%%%%%%%%%%%%% MAP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_x = [-27.75 50];
plot_y = [-36.8 7.9];

axes(axe1) % 2022
set(gca,'xaxislocation','top','FontSize',fontsize)
text(-29.75, 7, 'a','FontWeight','bold','FontSize',fontsize+3);

% ice-sheet bed contour map, colorbar
surf(BMv5_for_nevis_catchment.X_km/1e3,BMv5_for_nevis_catchment.Y_km./1e3,...
   BMv5_for_nevis_catchment.B_km,'EdgeColor','None','Facealpha',0.5); % plot surface below station points
view([0 90])
caxis([-600 200]); colormap(bone)

t2=colorbar('south');
hold all
set(get(t2,'xlabel'),'String',{'Bed Elevation  [ m a.s.l. ]'},'FontSize',fontsize+1,'Fontname','Avenir');
set(t2, 'Position', [0.05 0.6448 0.18 0.003],'XTick',-600:200:200,...
   'XTickLabel',[-600,-400,-200,0,200],'TickDirection','both');

hold all
% surface ice sheet contours
surf_contours = 0:100:1900; % m a.s.l.
[C, h] = contour(BMv5_for_nevis_catchment.X_km./1e3, BMv5_for_nevis_catchment.Y_km./1e3, ...
    BMv5_for_nevis_catchment.S_km, surf_contours,'Color',[0.7 0.7 0.7],'LineWidth',0.8);
clabel(C,h,'FontSize',7,'Color',[0.7 0.7 0.7],'FontName','Avenir','LabelSpacing',400)

% main array GNSS
plot3(xy_sta_22_short(:,1),xy_sta_22_short(:,2),5e3.*ones(length(xy_sta_22_short),1),'k^',...
    'MarkerSize',TriangleSize-2,'MarkerFaceColor','w','markerEdgeColor',metal,...
    'LineWidth',1.2);

% ROI that's relevant to the GNSS array = [-30 100] in X; [-40 10] in Y
% if the lake location is not within relevant catchment, throw it out !! 
ROI_relevant_xv = [-15 60 104 104 15 -10 -33 -33 -15]; % km in nevis-model space
ROI_relevant_yv = [15 4 -8 -50 -50 -41 -33 -1 15];
plot3(ROI_relevant_xv,ROI_relevant_yv,1e4.*[1 1 1 1 1 1 1 1 1],'-','LineWidth',2,'Color',[0.6 0.2 0]);

% lakes: FASTER S2 2022
[HF_lakes_to_plot_22] = find(environs_lakes_2022.laketypeing_dates(:,7)==1); % just HF lakes   
[MU_lakes_to_plot_22] = find(environs_lakes_2022.laketypeing_dates(:,7)==2); % just moulin-draining lakes 
[lakes_to_plot_22] = find(~isnan(environs_lakes_2022.laketypeing_dates(:,7))); % all good lakes   
[NE_lakes_to_plot_22] = find(environs_lakes_2022.laketypeing_dates(:,7)==5); % just no-exit freezers
[OS_lakes_to_plot_22] = find(environs_lakes_2022.laketypeing_dates(:,7)==3); % just overspillers
for i=1:length(OS_lakes_to_plot_22)
    IN_OS(i,1) = inpolygon(environs_lakes_2022.X_km(OS_lakes_to_plot_22(i)),environs_lakes_2022.Y_km(OS_lakes_to_plot_22(i)),...
        ROI_relevant_xv, ROI_relevant_yv);
end
[OS_lakes_to_plot_22_INPOLY] = OS_lakes_to_plot_22(IN_OS);

% plotting in order for legend
fill3(environs_lakes_2022.boundaries(lakes_to_plot_22(1)).XY_km_local(:,1), ...
        environs_lakes_2022.boundaries(lakes_to_plot_22(1)).XY_km_local(:,2), ...
        3000.*ones(length(environs_lakes_2022.boundaries(lakes_to_plot_22(1)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',lake_blue(1:3));
fill3(environs_lakes_2022.boundaries(HF_lakes_to_plot_22(1)).XY_km_local(:,1), ...
        environs_lakes_2022.boundaries(HF_lakes_to_plot_22(1)).XY_km_local(:,2), ...
        3000.*ones(length(environs_lakes_2022.boundaries(HF_lakes_to_plot_22(1)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',goldenrod);
fill3(environs_lakes_2022.boundaries(MU_lakes_to_plot_22(1)).XY_km_local(:,1), ...
        environs_lakes_2022.boundaries(MU_lakes_to_plot_22(1)).XY_km_local(:,2), ...
        3000.*ones(length(environs_lakes_2022.boundaries(MU_lakes_to_plot_22(1)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',[0.4 0.8 0.8]);
fill3(environs_lakes_2022.boundaries(NE_lakes_to_plot_22(1)).XY_km_local(:,1), ...
        environs_lakes_2022.boundaries(NE_lakes_to_plot_22(1)).XY_km_local(:,2), ...
        3000.*ones(length(environs_lakes_2022.boundaries(NE_lakes_to_plot_22(1)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',plum(1:3));

% Natalie Turner's lake connections for overspill lakes, JUST IN ROI
for i=20
    if strcmp(result{i,3}, 'L')
    plot3([environs_lakes_2022.X_km(i) environs_lakes_2022.X_km(result{i,4})],...
          [environs_lakes_2022.Y_km(i) environs_lakes_2022.Y_km(result{i,4})],[300 300],...
          '-','LineWidth',0.6,'Color',[0 0.55 1])
    else end
end

% cluster cOneXiOnS
grey_connections = lake_blue(1:3).*1.3; 
% 10s
cluster_10s = [34, 47, 72, 78, 81];
cc1 = plot3([environs_lakes_2022.X_km(cluster_10s(1)) environs_lakes_2022.X_km(cluster_10s(2))],...
    [environs_lakes_2022.Y_km(cluster_10s(1)) environs_lakes_2022.Y_km(cluster_10s(2))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
cc11 =plot3([environs_lakes_2022.X_km(cluster_10s(2)) environs_lakes_2022.X_km(cluster_10s(3))],...
    [environs_lakes_2022.Y_km(cluster_10s(2)) environs_lakes_2022.Y_km(cluster_10s(3))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
cc12 =plot3([environs_lakes_2022.X_km(cluster_10s(3)) environs_lakes_2022.X_km(cluster_10s(4))],...
    [environs_lakes_2022.Y_km(cluster_10s(3)) environs_lakes_2022.Y_km(cluster_10s(4))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
cc13 =plot3([environs_lakes_2022.X_km(cluster_10s(4)) environs_lakes_2022.X_km(cluster_10s(5))],...
    [environs_lakes_2022.Y_km(cluster_10s(4)) environs_lakes_2022.Y_km(cluster_10s(5))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
text(-17.5, -16,5e3, 'C1','FontSize',fontsize+7,'FontName','AvenirBold','FontAngle','italic',...
    'FontWeight','bold','Color',grey_connections);
% 20s
cluster_20s = [143, 155, 156, 187];
cc2 = plot3([environs_lakes_2022.X_km(cluster_20s(1)) environs_lakes_2022.X_km(cluster_20s(3))],...
    [environs_lakes_2022.Y_km(cluster_20s(1)) environs_lakes_2022.Y_km(cluster_20s(3))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
cc21 =plot3([environs_lakes_2022.X_km(cluster_20s(2)) environs_lakes_2022.X_km(cluster_20s(3))],...
    [environs_lakes_2022.Y_km(cluster_20s(2)) environs_lakes_2022.Y_km(cluster_20s(3))],[2e3 2e3],...
    ':','LineWidth',2,'Color',grey_connections);
cc22 =plot3([environs_lakes_2022.X_km(cluster_20s(1)) environs_lakes_2022.X_km(cluster_20s(4))],...
    [environs_lakes_2022.Y_km(cluster_20s(1)) environs_lakes_2022.Y_km(cluster_20s(4))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
cc23 =plot3([environs_lakes_2022.X_km(cluster_20s(4)) environs_lakes_2022.X_km(cluster_20s(3))],...
    [environs_lakes_2022.Y_km(cluster_20s(4)) environs_lakes_2022.Y_km(cluster_20s(3))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
text(22, -18,5e3, 'C3','FontSize',fontsize+7,'FontName','AvenirBold','FontAngle','italic',...
    'FontWeight','bold','Color',grey_connections);
% #mid
cluster_mid = [168, 173, 183, 185, 201, 223];
cc3 = plot3([environs_lakes_2022.X_km(cluster_mid(1)) environs_lakes_2022.X_km(cluster_mid(2))],...
    [environs_lakes_2022.Y_km(cluster_mid(1)) environs_lakes_2022.Y_km(cluster_mid(2))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
cc31 =plot3([environs_lakes_2022.X_km(cluster_mid(2)) environs_lakes_2022.X_km(cluster_mid(3))],...
    [environs_lakes_2022.Y_km(cluster_mid(2)) environs_lakes_2022.Y_km(cluster_mid(3))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
cc32 =plot3([environs_lakes_2022.X_km(cluster_mid(4)) environs_lakes_2022.X_km(cluster_mid(5))],...
    [environs_lakes_2022.Y_km(cluster_mid(4)) environs_lakes_2022.Y_km(cluster_mid(5))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
cc33 =plot3([environs_lakes_2022.X_km(cluster_mid(6)) environs_lakes_2022.X_km(cluster_mid(5))],...
    [environs_lakes_2022.Y_km(cluster_mid(6)) environs_lakes_2022.Y_km(cluster_mid(5))],[2e3 2e3],...
    '-','LineWidth',2,'Color',grey_connections);
text(24.5, -28.5,5e3, 'C2.S','FontSize',fontsize+7,'FontName','Avenirbold','FontAngle','italic',...
    'FontWeight','bold','Color',grey_connections);
text(33.5, 3.3,5e3, 'C2.N','FontSize',fontsize+7,'FontName','AvenirBold','FontAngle','italic',...
    'FontWeight','bold','Color',grey_connections);

% Natalie Turner's lake connections for overspill lakes, JUST IN ROI
for j=1:1:length(OS_lakes_to_plot_22_INPOLY)
    i = OS_lakes_to_plot_22_INPOLY(j);
    if strcmp(result{i,3}, 'L')
    plot3([environs_lakes_2022.X_km(i) environs_lakes_2022.X_km(result{i,4})],...
          [environs_lakes_2022.Y_km(i) environs_lakes_2022.Y_km(result{i,4})],[300 300],...
          '-','LineWidth',0.4,'Color',[0 0.55 1])
    else end
end
for j=1:1:length(OS_lakes_to_plot_22_INPOLY)
    i = OS_lakes_to_plot_22_INPOLY(j);
if strcmp(result{i,3}, 'M')
    plot3([environs_lakes_2022.X_km(i) moulins_x_2022(result{i,4})],...
        [environs_lakes_2022.Y_km(i) moulins_y_2022(result{i,4})],[300 300],...
        '-','LineWidth',0.4,'Color',[0 0.55 1])
    % plot the moulin
    plot3(moulins_x_2022(result{i,4}), moulins_y_2022(result{i,4}),...
        300.*ones(length(moulins_x_2022(result{i,4})),1),'ko',...
        'MarkerSize',TriangleSize-4,'MarkerFaceColor','w','markerEdgeColor','k',...
        'LineWidth',0.5);
    else end
end

% fill those lakes, girl
for i=1:1:length(lakes_to_plot_22)
    fill3(environs_lakes_2022.boundaries(lakes_to_plot_22(i)).XY_km_local(:,1), ...
        environs_lakes_2022.boundaries(lakes_to_plot_22(i)).XY_km_local(:,2), ...
        3000.*ones(length(environs_lakes_2022.boundaries(lakes_to_plot_22(i)).XY_km_local(:,2)),1),... 
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',lake_blue(1:3));
end
for i=1:1:length(HF_lakes_to_plot_22)
    fill3(environs_lakes_2022.boundaries(HF_lakes_to_plot_22(i)).XY_km_local(:,1), ...
        environs_lakes_2022.boundaries(HF_lakes_to_plot_22(i)).XY_km_local(:,2), ...
        3000.*ones(length(environs_lakes_2022.boundaries(HF_lakes_to_plot_22(i)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',goldenrod);
end
for i=1:1:length(MU_lakes_to_plot_22)
    fill3(environs_lakes_2022.boundaries(MU_lakes_to_plot_22(i)).XY_km_local(:,1), ...
        environs_lakes_2022.boundaries(MU_lakes_to_plot_22(i)).XY_km_local(:,2), ...
         3000.*ones(length(environs_lakes_2022.boundaries(MU_lakes_to_plot_22(i)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',[0.4 0.8 0.8]);
end
for i=1:1:length(NE_lakes_to_plot_22)
    fill3(environs_lakes_2022.boundaries(NE_lakes_to_plot_22(i)).XY_km_local(:,1), ...
        environs_lakes_2022.boundaries(NE_lakes_to_plot_22(i)).XY_km_local(:,2), ...
        3000.*ones(length(environs_lakes_2022.boundaries(NE_lakes_to_plot_22(i)).XY_km_local(:,2)),1),...
        'o','MarkerSize',1,'EdgeColor','none','FaceColor',plum(1:3));
end

% HF dates of drainages range (plot as text)
environs_lakes_2022.laketypeing_dates(34,3:4) = [194, 195];
for i=1:1:length(HF_lakes_to_plot_22)
  j=HF_lakes_to_plot_22(i);
    if environs_lakes_2022.laketypeing_dates(j,3)+1 == environs_lakes_2022.laketypeing_dates(j,4)
      t=text(environs_lakes_2022.X_km(j)-1.0, environs_lakes_2022.Y_km(j),5e3,...
        sprintf('%d',environs_lakes_2022.laketypeing_dates(j,3)+1), ...
        'FontSize',fontsize-2,'FontName','Helvetica','FontWeight','bold');
      t.Clipping = "on";
   else
      t=text(environs_lakes_2022.X_km(j)-1.5, environs_lakes_2022.Y_km(j),5e3,...
        sprintf('%d\x2013%d',environs_lakes_2022.laketypeing_dates(j,3)+1,environs_lakes_2022.laketypeing_dates(j,4)), ...
        'FontSize',fontsize-2,'FontName','Helvetica','FontWeight','bold');
      t.Clipping = "on";
   end
end

% Moulin dates of drainages (plot as text)
for i=1:1:length(MU_lakes_to_plot_22)
   j=MU_lakes_to_plot_22(i);
   if environs_lakes_2022.laketypeing_dates(j,3)+1 == environs_lakes_2022.laketypeing_dates(j,4)
      t=text(environs_lakes_2022.X_km(j)-1.0, environs_lakes_2022.Y_km(j),5e3,...
        sprintf('%d',environs_lakes_2022.laketypeing_dates(j,3)+1), ...
        'FontSize',fontsize-2,'FontName','Helvetica','FontWeight','bold');
      t.Clipping = "on";
   else 
    t=text(environs_lakes_2022.X_km(j)-1.5, environs_lakes_2022.Y_km(j),5e3,...
        sprintf('%d\x2013%d',environs_lakes_2022.laketypeing_dates(j,3)+1,environs_lakes_2022.laketypeing_dates(j,4)), ...
        'FontSize',fontsize-2,'FontName','Helvetica','FontWeight','bold');
    t.Clipping = "on";
   end
end

% GNSS site names
for i=1; text(xy_sta_22_short(i,1)-3.75,xy_sta_22_short(i,2)+0.55,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=2; text(xy_sta_22_short(i,1)-4.25,xy_sta_22_short(i,2)+0.55,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=3; text(xy_sta_22_short(i,1)+0.75,xy_sta_22_short(i,2)+0.55,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=5; text(xy_sta_22_short(i,1)+1.25,xy_sta_22_short(i,2)-0.55,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=6; text(xy_sta_22_short(i,1)+1.00,xy_sta_22_short(i,2)+0.10,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=9; text(xy_sta_22_short(i,1)+0.20,xy_sta_22_short(i,2)+2.00,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=4; text(xy_sta_22_short(i,1)-2.85,xy_sta_22_short(i,2)-1.55,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic','Rotation',30); end
for i=7; text(xy_sta_22_short(i,1)-1.5,xy_sta_22_short(i,2)-1.25,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=8; text(xy_sta_22_short(i,1)-3.25,xy_sta_22_short(i,2)+0.75,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=10:11; text(xy_sta_22_short(i,1)-3.75,xy_sta_22_short(i,2)-0.1,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=12:13; text(xy_sta_22_short(i,1)+1.25,xy_sta_22_short(i,2)-0.55,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=14; text(xy_sta_22_short(i,1)-3.75,xy_sta_22_short(i,2)+0.55,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=15:16; text(xy_sta_22_short(i,1)+1.25,xy_sta_22_short(i,2)+0.55,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=17; text(xy_sta_22_short(i,1)-1.5,xy_sta_22_short(i,2)-1.0,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=18; text(xy_sta_22_short(i,1)-3.5,xy_sta_22_short(i,2)+0.25,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=19; text(xy_sta_22_short(i,1)-3.5,xy_sta_22_short(i,2)-0.5,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=20; text(xy_sta_22_short(i,1)+0.75,xy_sta_22_short(i,2)+0.55,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic'); end
for i=21; text(xy_sta_22_short(i,1)+0.50,xy_sta_22_short(i,2)-0.55,5e3,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-1,'FontAngle','italic','Rotation',-45); end
% lake names
text(2, 4 ,5e3, 'L1C','Fontname','AvenirBold','FontSize',fontsize+2,'FontAngle','italic','Color','k')
text(3.25, 0.75 ,5e3, 'L1A','Fontname','AvenirBold','FontSize',fontsize+2,'FontAngle','italic','Color','k')
text(1.25, -3, 5e3, 'L1B','Fontname','AvenirBold','FontSize',fontsize+2,'FontAngle','italic','Color','k')
text(24, -11.5, 5e3, 'L2A','Fontname','AvenirBold','FontSize',fontsize+2,'FontAngle','italic','Color','k')
text(13.5, -20 ,5e3, 'L2B','Fontname','AvenirBold','FontSize',fontsize+2,'FontAngle','italic','Color','k')
text(27, -17.5 ,5e3, 'L2C','Fontname','AvenirBold','FontSize',fontsize+2,'FontAngle','italic','Color','k')
text(38, -21.25, 5e3, 'L3A','Fontname','AvenirBold','FontSize',fontsize+2,'FontAngle','italic','Color','k')
text(31.5, -27 ,5e3, 'L3B','Fontname','AvenirBold','FontSize',fontsize+2,'FontAngle','italic','Color','k')

% scalebar
plot3([-10 0],[-35 -35],[5e3 5e3],'k-','LineWidth',2);
text(-7, -33.75, 5e3, '10 km','FontName','AvenirBold','FontSize',fontsize+1)

lgd = legend('Bed Elevation','Ice Elevation','GNSS','ROI',...
    'Overspill drainage','HF drainage','Moulin drainage','No-exit, frozen',...
    'Overspill connections','HF Clusters ({\itp}-value<0.05)');
set(lgd,'Position',[0.510 0.5825 0.01 0.01],'FontSize',fontsize,'FontName','Avenir',...
    'NumColumns',5,'Box','off','TextColor',lake_blue(1:3))

%axis equal; 
grid on; box on;
text(-26, 3.5, 5e3, '2022','FontSize',fontsize+4,'FontName','Helvetica','FontWeight','bold','Rotation',45)
set(gca,'FontName','Avenir','FontSize',fontsize,'xtick',-80:10:110,...
    'ytick',-50:10:30,'LineWidth',0.6)
xlim([plot_x(1) plot_x(2)]); ylim([plot_y(1) plot_y(2)]);
axe1.Layer = 'top';   axe1.ClippingStyle = "rectangle";  

%% print figure
print(gcf,'-dpng','-r300',sprintf('../../paperfigs/paperfig5_distillations_2022_260129_repository.png')); 