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
% LAS 2024-10-20: days of C6 cluster
% LAS 25-04-24: FASTER lake locations and 30-min velocities
% LAS 26-01-19: +error envelopes on strain rates and c_dot; repository

close all; clear all

% load runoff data
load('../../nevis_lakes_cluster/racmo_station_2023_300m_index.mat') % runoff_2023_nevis(:,ID(6,1))./10 = cm w.e. at SQ13
load('../../nevis_lakes_cluster/runoff_2023_nevis_noSK_300m.mat');
runoff_2023_nevis = runoff_2023_nevis.runoff;
% load station names
load('station_names_2023.mat');
num_sta = length(station_names); % number of stations

% load c_dot timeseries -- loose
load('daily_epsilon_zz_2023R_30min_BF2_UP4_sZERO_OUT6_clean_w36_t12_260119.mat') % 30 min vels; No smooth; window=36h; threshold=8h
daily_epsilon_zz_loose = daily_epsilon_zz;

% load c_dot timeseries -- tight
load('daily_epsilon_zz_2023R_30min_BF2_UP4_sZERO_OUT6_clean_w18_t6_260119.mat')  % 30-min vels; No smooth; window=18h; threshold=4h
daily_epsilon_zz_tight = daily_epsilon_zz;

% load strain rates -- loose
load('daily_strain_rates_2023R_30min_BF2_UP4_sZERO_OUT6_clean_w36_t12_260119.mat') % 30-min vels; No smooth; window=36h; threshold=8h
daily_strain_rates_loose = daily_strain_rates_2023;
strain_time = daily_strain_rates_2023.time;
station_names_short = upper(daily_strain_rates_2023.station_names); % less SQ33

% load strain rates -- tight
load('daily_strain_rates_2023R_30min_BF2_UP4_sZERO_OUT6_clean_w18_t6_260119.mat') % 30-min vels; No smooth; window=18h; threshold=4h
daily_strain_rates_tight = daily_strain_rates_2023; 

% Tight constraints at station-days 2023:
tight_100s = 188:196;
tight_200s = 195:205;
tight_300s = 195:205;
tight_MLOW = 188:205;

% % tight constraints on all days
% tight_100s = 150:250;
% tight_200s = 150:250;
% tight_300s = 150:250;
% tight_MLOW = 150:250;
% tight_QIET = 150:250;
% tight_MHIH = 150:250;

%% COMBINE TIGHT AND LOOSE RECORDS FOR CDOT
% make timeseries logical for where to use 'tight' strain rates

for i=1:3 % MLOW, MHIH, QIET
    % first fill with loose constraints (longer sliding LS window)
    % loose daily_epsilon_zz
    daily_epsilon_zz(i).t22 = daily_epsilon_zz_loose(i).t22; % time
    daily_epsilon_zz(i).c_dot_delta_t_combo = daily_epsilon_zz_loose(i).c_dot_delta_t; % c_dot
    daily_epsilon_zz(i).c_dot_delta_t_err_combo = daily_epsilon_zz_loose(i).c_dot_delta_t_err; % c_dot_err
    for j=1:1:length(daily_epsilon_zz(i).t22)
        if daily_epsilon_zz_tight(i).t22(j,1) >= tight_MLOW(1) && daily_epsilon_zz_tight(i).t22(j,1) <= tight_MLOW(end)+1
        % assign outputs from tight constraint (shorter sliding LS window)
        daily_epsilon_zz(i).c_dot_delta_t_combo(j,1) = daily_epsilon_zz_tight(i).c_dot_delta_t(j,1); % c_dot tight
        daily_epsilon_zz(i).c_dot_delta_t_combo_err(j,1) = daily_epsilon_zz_tight(i).c_dot_delta_t_err(j,1); % c_dot_err tight
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
    daily_epsilon_zz(i).c_dot_delta_t_err_combo = daily_epsilon_zz_loose(i).c_dot_delta_t_err; % c_dot_err
    for j=1:1:length(daily_epsilon_zz(i).t22)
        if daily_epsilon_zz_tight(i).t22(j,1) >= tight_100s(1) && daily_epsilon_zz_tight(i).t22(j,1) <= tight_100s(end)+1
        % assign outputs from tight constraint
        daily_epsilon_zz(i).c_dot_delta_t_combo(j,1) = daily_epsilon_zz_tight(i).c_dot_delta_t(j,1); % c_dot tight
        daily_epsilon_zz(i).c_dot_delta_t_combo_err(j,1) = daily_epsilon_zz_tight(i).c_dot_delta_t_err(j,1); % c_dot_err tight
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
    daily_epsilon_zz(i).c_dot_delta_t_err_combo = daily_epsilon_zz_loose(i).c_dot_delta_t_err; % c_dot_err
    for j=1:1:length(daily_epsilon_zz(i).t22)
        if daily_epsilon_zz_tight(i).t22(j,1) >= tight_200s(1) && daily_epsilon_zz_tight(i).t22(j,1) <= tight_200s(end)+1
        % assign outputs from tight constraint
        daily_epsilon_zz(i).c_dot_delta_t_combo(j,1) = daily_epsilon_zz_tight(i).c_dot_delta_t(j,1); % c_dot tight
        daily_epsilon_zz(i).c_dot_delta_t_combo_err(j,1) = daily_epsilon_zz_tight(i).c_dot_delta_t_err(j,1); % c_dot_err tight
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
    daily_epsilon_zz(i).c_dot_delta_t_err_combo = daily_epsilon_zz_loose(i).c_dot_delta_t_err; % c_dot_err
    for j=1:1:length(daily_epsilon_zz(i).t22)
        if daily_epsilon_zz_tight(i).t22(j,1) >= tight_300s(1) && daily_epsilon_zz_tight(i).t22(j,1) <= tight_300s(end)+1
        % assign outputs from tight constraint
        daily_epsilon_zz(i).c_dot_delta_t_combo(j,1) = daily_epsilon_zz_tight(i).c_dot_delta_t(j,1); % c_dot tight
        daily_epsilon_zz(i).c_dot_delta_t_combo_err(j,1) = daily_epsilon_zz_tight(i).c_dot_delta_t_err(j,1); % c_dot_err tight
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


%% plotting preliminaries 
% 2023 lake dates
lake_dates = [191.2, 196.75, 202.85, 200.50, 195.25]; % 100s, MHIH (IDs 150&152), 200s, 300s, IDs 154&162
lake_fill_dates = [179, 190, 183, 190, 190]; % 100s, MHIH, 200s, 300s, IDs 154&162
% strain-rate time
delta_t = daily_epsilon_zz(1).delta_t;
% racmo
racmo_time = 1.5:1:365.5; 
% t_0 for strain_time
[ID_t0] = find(strain_time>=165.01, 1, 'first');
[ID_165] = find(strain_time>=165.01, 1, 'first');
[ID_198] = find(strain_time>=198.0, 1, 'first');
[ID_200] = find(strain_time>=200.0, 1, 'first');

%% figure -- 2023 stack: c_dot, strain rate, runoff
figure(1); clf; 
% 5-panel stack:
set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',1.15.*[1 0 18.3 24.7]);
axe1 = axes('Position',[0.08 0.811 0.855 0.182],'Box','on','NextPlot','add','XTickLabels',[]);
axe2 = axes('Position',[0.08 0.621 0.855 0.182],'Box','on','NextPlot','add','XTickLabels',[]);
axe3 = axes('Position',[0.08 0.426 0.855 0.182],'Box','on','NextPlot','add','XTickLabels',[]);
axe4 = axes('Position',[0.08 0.231 0.855 0.182],'Box','on','NextPlot','add','XTickLabels',[]);
axe5 = axes('Position',[0.08 0.036 0.855 0.182],'Box','on','NextPlot','add');

t_plot = [197 206];

y_plot_strain = [-0.05 0.10]; % 2023 cluster C4
y_plot_strain_tick = 0.025;
y_plot_cdot = [-3.75 0.75];
y_plot_cdot_tick = 0.75; y_plot_cdot_interval = y_plot_cdot_tick;
plot_cdot_range = 6; %y_plot_cdot(2)-y_plot_cdot(1);
run_c_dot_label = [{'2'};{'4'};{'6'};{'-0.75'};{'0.0'};{'0.75'}];

fontsize = 9; 
TriangleSize = 7;
DiamondSize = 1.5;
plum = [211, 160, 211, 100]./255;
lake_blue = [0    0.4470    0.7410   0.06]; % fourth value is FaceAlpha built in!
% dock at eel pond 
sky_blue = [111, 169, 228]./255; % SNL
metal = [87, 115, 131]./255; % NNL
oar = [251, 219, 154]./255; % NENL
handle = [161, 37, 49]./255; % NL
dark_oar = [164, 114, 63]./255;

% % ORDER IN TRACK READ-IN:
% % 1 = MHIH; 2 = MLOW; 3=QIET; 
% % 4=SQ11;   5=SQ12;   6=SQ13;   7=SQ14;   8=SQ15;   9=SQ16;
% % 10=SQ21;  11=SQ22;  12=SQ23;  13=SQ24;  14=SQ25;  15=SQ26;  16=SQ27;
% % 17=SQ31;  18=SQ32;  19=SQ34;  20=SQ35;  21=SQ36;  22=SQ37;

%%%%%%%%%%%%%%%%%%%%%%%%% 100s %%%%%%%%%%%%%%%%%%%%%%%%% 
axes(axe1) 
colororder({'k','k'})
yyaxis right
% \dot{c}
for i=4:1:9
    % c_dot_err
    tt = [daily_epsilon_zz(i).t22; flipud(daily_epsilon_zz(i).t22)];
    ee = [daily_epsilon_zz(i).c_dot_delta_t_combo+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165))+(3.*daily_epsilon_zz(i).c_dot_delta_t_err_combo);...
        flipud(daily_epsilon_zz(i).c_dot_delta_t_combo+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165))-(3.*daily_epsilon_zz(i).c_dot_delta_t_err_combo))];
    nanIdx = isnan(ee); ee(nanIdx) = 0;
    patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
end
% \dot{c}
for i=4:1:9
    plot(daily_epsilon_zz(i).t22, ...
        (daily_epsilon_zz(i).c_dot_delta_t_combo)+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165)),...
        '-','LineWidth',1.5); hold on;
end
colororder('default');  CMap = get(gca, 'ColorOrder');
% runoff
patch_vert = ((4*y_plot_cdot_interval)/plot_cdot_range).*(nanmean(runoff_2023_nevis(:,ID(4:9,1)),2)./10); % shrink to c_dot scale [(num_intervals*interval_size)/RUnoff Range (8 cm)]
% patch_x1 = racmo_time-0.5+0.25; patch_x2 = patch_x1+0.5;
patch_x1 = racmo_time-0.5; patch_x2 = patch_x1+1;
patch_y1 = ones(length(racmo_time),1).*y_plot_cdot(1);
patch_y2 = patch_y1 + patch_vert; 
for i=1:1:length(racmo_time)
    patch([patch_x1(i) patch_x2(i) patch_x2(i) patch_x1(i)],...
          [patch_y1(i) patch_y1(i) patch_y2(i) patch_y2(i)],...
          [0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.20); hold on;
end
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_cdot(1) y_plot_cdot(2)])
ylabel('Runoff [cm w.e.] \hspace{6mm} $\dot{c}$  [ m d$^{-1}$ ] \hspace{0mm}','Interpreter','latex','FontSize',fontsize)
set(gca,'ytick',y_plot_cdot(1)+y_plot_cdot_tick:y_plot_cdot_tick:y_plot_cdot(2),'YTickLabel',run_c_dot_label);
text(t_plot(2)-0.15,y_plot_cdot(2)-1.10,'within-basin 950s','FontName','AvenirBold','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','right')
text(t_plot(1)-0.80, y_plot_cdot(2),'a','FontName','Helvetica','FontSize',fontsize+2);

text(t_plot(1)+0.5, y_plot_cdot(2)-0.20, 'SQ11','FontName','AvenirBold','FontSize',fontsize,'Color',[0    0.4470    0.7410],'HorizontalAlignment','center');
text(t_plot(1)+1.5, y_plot_cdot(2)-0.20, 'SQ12','FontName','AvenirBold','FontSize',fontsize,'Color',[0.8500    0.3250    0.0980],'HorizontalAlignment','center');
text(t_plot(1)+2.5, y_plot_cdot(2)-0.20, 'SQ13','FontName','AvenirBold','FontSize',fontsize,'Color',[0.9290    0.6940    0.1250],'HorizontalAlignment','center');
text(t_plot(1)+3.5, y_plot_cdot(2)-0.20, 'SQ14','FontName','AvenirBold','FontSize',fontsize,'Color',[0.4940    0.1840    0.5560],'HorizontalAlignment','center');
text(t_plot(1)+4.5, y_plot_cdot(2)-0.20, 'SQ15','FontName','AvenirBold','FontSize',fontsize,'Color',[0.4660    0.6740    0.1880],'HorizontalAlignment','center');
text(t_plot(1)+5.5, y_plot_cdot(2)-0.20, 'SQ16','FontName','AvenirBold','FontSize',fontsize,'Color',[0.3010    0.7450    0.9330],'HorizontalAlignment','center');

text(lake_dates(1)+0.25, y_plot_cdot(1)+0.20, {'L1A';'L1B'},'FontName','AvenirBold','FontSize',fontsize+1,...
    'FontAngle','italic','HorizontalAlignment','left','Rotation',0); 

yyaxis left % \dot{epsilon}_{lon}
for i=1 % lakes filling and draining
    rectangle('Position',[lake_fill_dates(i) y_plot_strain(1) lake_dates(i)-lake_fill_dates(i) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',lake_blue,'EdgeColor','none'); 
    rectangle('Position',[lake_dates(i) y_plot_strain(1) 1 y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[1 1 0.7],'EdgeColor','none');
end
for i=2:1:length(lake_dates)-1 % other drainages
    rectangle('Position',[lake_dates(i) y_plot_strain(1) 1 y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[0.90 0.90 0.90],'EdgeColor','none');
end
for i=length(lake_dates) % 195 drainage
    rectangle('Position',[lake_dates(i) y_plot_strain(1) 1 y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[1 1 0.7],'EdgeColor','none');
end

% strain rates
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(4,5,:)-daily_strain_rates_combo.lon_yr_100s(4,5,ID_t0)),'d-','MarkerSize',DiamondSize,'Color',CMap(1,:)); hold on; %SQ11-12
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(6,7,:)-daily_strain_rates_combo.lon_yr_100s(6,7,ID_t0)),'d-','MarkerSize',DiamondSize,'Color',CMap(2,:)); %SQ13-14
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(8,9,:)-daily_strain_rates_combo.lon_yr_100s(8,9,ID_t0)),'d-','MarkerSize',DiamondSize,'Color',CMap(3,:)); %SQ15-16
% strain-rate error envelopes (any data gaps --> error estimate == 0)
tt = [strain_time; flipud(strain_time)];
ee = [squeeze(daily_strain_rates_combo.lon_yr_100s(4,5,:)+3.*daily_strain_rates_combo.delta_lon_yr_100s(4,5,:)-daily_strain_rates_combo.lon_yr_100s(4,5,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_100s(4,5,:)-3.*daily_strain_rates_combo.delta_lon_yr_100s(4,5,:)-daily_strain_rates_combo.lon_yr_100s(4,5,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_100s(6,7,:)+3.*daily_strain_rates_combo.delta_lon_yr_100s(6,7,:)-daily_strain_rates_combo.lon_yr_100s(6,7,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_100s(6,7,:)-3.*daily_strain_rates_combo.delta_lon_yr_100s(6,7,:)-daily_strain_rates_combo.lon_yr_100s(6,7,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_100s(8,9,:)+3.*daily_strain_rates_combo.delta_lon_yr_100s(8,9,:)-daily_strain_rates_combo.lon_yr_100s(8,9,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_100s(8,9,:)-3.*daily_strain_rates_combo.delta_lon_yr_100s(8,9,:)-daily_strain_rates_combo.lon_yr_100s(8,9,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5); hold on;
% strain rates
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(4,5,:)-daily_strain_rates_combo.lon_yr_100s(4,5,ID_t0)),'d-','MarkerSize',DiamondSize,'Color',CMap(1,:)); hold on; %SQ11-12
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(6,7,:)-daily_strain_rates_combo.lon_yr_100s(6,7,ID_t0)),'d-','MarkerSize',DiamondSize,'Color',CMap(2,:)); %SQ13-14
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(8,9,:)-daily_strain_rates_combo.lon_yr_100s(8,9,ID_t0)),'d-','MarkerSize',DiamondSize,'Color',CMap(3,:)); %SQ15-16

colororder('default');
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_strain(1) y_plot_strain(2)])
ylabel('$\dot{\epsilon}_{lon}(t)-\dot{\epsilon}_{lon}(t_{165})$  [ yr$^{-1}$ ]','Interpreter','latex','FontSize',fontsize+1)
legend(['SQ11' char(8211) '12'],['SQ13' char(8211) '14'],['SQ15' char(8211) '16'],'NumColumns',3,...
    'Location','West','EdgeColor','none','color','none')
    % 'Position',[0.205 0.805 0.01 0.01],'EdgeColor','none')
set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.6,'FontSize',fontsize); 
set(gca,'ytick',y_plot_strain(1)+y_plot_strain_tick:y_plot_strain_tick:y_plot_strain(2))
grid on; box on; axe1.Layer = 'top';


%%%%%%%%%%%%%%%%%%%%%%%%% 100s-to-Ties %%%%%%%%%%%%%%%%%%%%%%%%% 
axes(axe2) 
colororder({'k','k'})
yyaxis right 
for i=2:3
    % c_dot_err
    tt = [daily_epsilon_zz(i).t22; flipud(daily_epsilon_zz(i).t22)];
    ee = [daily_epsilon_zz(i).c_dot_delta_t_combo+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165))+(3.*daily_epsilon_zz(i).c_dot_delta_t_err_combo);...
        flipud(daily_epsilon_zz(i).c_dot_delta_t_combo+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165))-(3.*daily_epsilon_zz(i).c_dot_delta_t_err_combo))];
    nanIdx = isnan(ee); ee(nanIdx) = 0;
    patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
end
% \dot{c}
plot(daily_epsilon_zz(3).t22, (daily_epsilon_zz(3).c_dot_delta_t_combo)+(-1.*daily_epsilon_zz(3).c_dot_delta_t_combo(ID_165)),'-','LineWidth',1.5); hold on;
plot(daily_epsilon_zz(2).t22, (daily_epsilon_zz(2).c_dot_delta_t_combo)+(-1.*daily_epsilon_zz(2).c_dot_delta_t_combo(ID_165)),'-','LineWidth',1.5); hold on;
colororder('default');
% runoff
patch_vert = ((4*y_plot_cdot_interval)/plot_cdot_range).*(nanmean(runoff_2023_nevis(:,ID(2:3,1)),2)./10); % shrink to c_dot scale
% patch_x1 = racmo_time-0.5+0.25; patch_x2 = patch_x1+0.5;
patch_x1 = racmo_time-0.5; patch_x2 = patch_x1+1;
patch_y1 = ones(length(racmo_time),1).*y_plot_cdot(1);
patch_y2 = patch_y1 + patch_vert; 
for i=1:1:length(racmo_time)
    patch([patch_x1(i) patch_x2(i) patch_x2(i) patch_x1(i)],...
          [patch_y1(i) patch_y1(i) patch_y2(i) patch_y2(i)],...
          [0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.20); hold on;
end
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_cdot(1) y_plot_cdot(2)])
ylabel('Runoff [cm w.e.] \hspace{6mm} $\dot{c}$  [ m d$^{-1}$ ] \hspace{0mm}','Interpreter','latex','FontSize',fontsize)
set(gca,'ytick',y_plot_cdot(1)+y_plot_cdot_tick:y_plot_cdot_tick:y_plot_cdot(2),'YTickLabel',run_c_dot_label);
text(t_plot(2)-0.15,y_plot_cdot(2)-1.10,'950s-to-tiepoints','FontName','AvenirBold','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','right')
text(t_plot(1)-0.80, y_plot_cdot(2),'b','FontName','Helvetica','FontSize',fontsize+2);
text(t_plot(1)+0.5, y_plot_cdot(2)-0.20, 'QIET','FontName','AvenirBold','FontSize',fontsize,'Color',[0    0.4470    0.7410],'HorizontalAlignment','center');
text(t_plot(1)+1.5, y_plot_cdot(2)-0.20, 'MLOW','FontName','AvenirBold','FontSize',fontsize,'Color',[0.8500    0.3250    0.0980],'HorizontalAlignment','center');


yyaxis left % \dot{epsilon}_{lon}

for i=1:1:length(lake_dates) % other drainages
    rectangle('Position',[lake_dates(i) y_plot_strain(1) 1 y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[0.90 0.90 0.90],'EdgeColor','none');
end

% Ties
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(3,5,:)-daily_strain_rates_combo.lon_yr_100s(3,5,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(1,:)); hold on; %SQ12-QIET
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(3,9,:)-daily_strain_rates_combo.lon_yr_100s(3,9,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(2,:)); %SQ16-QIET
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(2,5,:)-daily_strain_rates_combo.lon_yr_100s(2,5,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(3,:)); %SQ12-MLOW
% strain-rate error envelopes (any data gaps --> error estimate == 0)
tt = [strain_time; flipud(strain_time)];
ee = [squeeze(daily_strain_rates_combo.lon_yr_100s(3,5,:)+3.*daily_strain_rates_combo.delta_lon_yr_100s(3,5,:)-daily_strain_rates_combo.lon_yr_100s(3,5,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_100s(3,5,:)-3.*daily_strain_rates_combo.delta_lon_yr_100s(3,5,:)-daily_strain_rates_combo.lon_yr_100s(3,5,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_100s(3,9,:)+3.*daily_strain_rates_combo.delta_lon_yr_100s(3,9,:)-daily_strain_rates_combo.lon_yr_100s(3,9,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_100s(3,9,:)-3.*daily_strain_rates_combo.delta_lon_yr_100s(3,9,:)-daily_strain_rates_combo.lon_yr_100s(3,9,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_100s(2,5,:)+3.*daily_strain_rates_combo.delta_lon_yr_100s(2,5,:)-daily_strain_rates_combo.lon_yr_100s(2,5,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_100s(2,5,:)-3.*daily_strain_rates_combo.delta_lon_yr_100s(2,5,:)-daily_strain_rates_combo.lon_yr_100s(2,5,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5); hold on;
% Ties
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(3,5,:)-daily_strain_rates_combo.lon_yr_100s(3,5,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(1,:)); hold on; %SQ12-QIET
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(3,9,:)-daily_strain_rates_combo.lon_yr_100s(3,9,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(2,:)); %SQ16-QIET
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(2,5,:)-daily_strain_rates_combo.lon_yr_100s(2,5,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(3,:)); %SQ12-MLOW

colororder('default');
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_strain(1) y_plot_strain(2)])
ylabel('$\dot{\epsilon}_{lon}(t)-\dot{\epsilon}_{lon}(t_{165})$  [ yr$^{-1}$ ]','Interpreter','latex','FontSize',fontsize+1)
legend(['SQ13' char(8211) 'QIET'],['SQ16' char(8211) 'QIET'],['SQ12' char(8211) 'MLOW'],'NumColumns',3,...
    'Location','West','EdgeColor','none','color','none');    % 'Position', [0.2275 0.630 0.01 0.01],'EdgeColor','none')
set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.6,'FontSize',fontsize); 
set(gca,'ytick',y_plot_strain(1)+y_plot_strain_tick:y_plot_strain_tick:y_plot_strain(2))
grid on; box on; axe2.Layer = 'top';

%%%%%%%%%%%%%%%%%%%%%%%%% Ties-to-200s %%%%%%%%%%%%%%%%%%%%%%%%% 
axes(axe3) 
colororder({'k','k'})
yyaxis right
for i=1:3
    % c_dot_err
    tt = [daily_epsilon_zz(i).t22; flipud(daily_epsilon_zz(i).t22)];
    ee = [daily_epsilon_zz(i).c_dot_delta_t_combo+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165))+(3.*daily_epsilon_zz(i).c_dot_delta_t_err_combo);...
        flipud(daily_epsilon_zz(i).c_dot_delta_t_combo+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165))-(3.*daily_epsilon_zz(i).c_dot_delta_t_err_combo))];
    nanIdx = isnan(ee); ee(nanIdx) = 0;
    patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
end
% \dot{c}
plot(daily_epsilon_zz(3).t22, (daily_epsilon_zz(3).c_dot_delta_t_combo)+(-1.*daily_epsilon_zz(3).c_dot_delta_t_combo(ID_165)),'-','LineWidth',1.5); hold on;
plot(daily_epsilon_zz(2).t22, (daily_epsilon_zz(2).c_dot_delta_t_combo)+(-1.*daily_epsilon_zz(2).c_dot_delta_t_combo(ID_165)),'-','LineWidth',1.5); hold on;
plot(daily_epsilon_zz(1).t22, (daily_epsilon_zz(1).c_dot_delta_t_combo)+(-1.*daily_epsilon_zz(1).c_dot_delta_t_combo(ID_165)),'-','LineWidth',1.5); hold on;
colororder('default');
% runoff
patch_vert = ((4*y_plot_cdot_interval)/plot_cdot_range).*(nanmean(runoff_2023_nevis(:,ID(2:3,1)),2)./10); % shrink to c_dot scale
% patch_x1 = racmo_time-0.5+0.25; patch_x2 = patch_x1+0.5;
patch_x1 = racmo_time-0.5; patch_x2 = patch_x1+1;
patch_y1 = ones(length(racmo_time),1).*y_plot_cdot(1);
patch_y2 = patch_y1 + patch_vert; 
for i=1:1:length(racmo_time)
    patch([patch_x1(i) patch_x2(i) patch_x2(i) patch_x1(i)],...
          [patch_y1(i) patch_y1(i) patch_y2(i) patch_y2(i)],...
          [0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.20); hold on;
end
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_cdot(1) y_plot_cdot(2)])
ylabel('Runoff [cm w.e.] \hspace{4mm} $\dot{c}$  [ m d$^{-1}$ ] \hspace{2mm}','Interpreter','latex','FontSize',fontsize)
set(gca,'ytick',y_plot_cdot(1)+y_plot_cdot_tick:y_plot_cdot_tick:y_plot_cdot(2),'YTickLabel',run_c_dot_label);
text(t_plot(2)-0.15,y_plot_cdot(2)-1.10,'tiepoints-to-1150s','FontName','AvenirBold','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','right')
text(t_plot(1)-0.80, y_plot_cdot(2),'c','FontName','Helvetica','FontSize',fontsize+2);
text(t_plot(1)+0.5, y_plot_cdot(2)-0.20, 'QIET','FontName','AvenirBold','FontSize',fontsize,'Color',[0    0.4470    0.7410],'HorizontalAlignment','center');
text(t_plot(1)+1.5, y_plot_cdot(2)-0.20, 'MLOW','FontName','AvenirBold','FontSize',fontsize,'Color',[0.8500    0.3250    0.0980],'HorizontalAlignment','center');
text(t_plot(1)+2.5, y_plot_cdot(2)-0.20, 'MHIH','FontName','AvenirBold','FontSize',fontsize,'Color',[0.9290    0.6940    0.1250],'HorizontalAlignment','center');

yyaxis left % \dot{epsilon}_{lon}

for i=1:1:length(lake_dates) % other drainages
    rectangle('Position',[lake_dates(i) y_plot_strain(1) 1 y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[0.90 0.90 0.90],'EdgeColor','none');
end
for i=2 % other drainages
    rectangle('Position',[lake_dates(i) y_plot_strain(1) 1 y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[1 1 0.7],'EdgeColor','none');
end
text(lake_dates(2)+0.30, y_plot_strain(1)+0.014, {'C5';'Southern'},'FontName','AvenirBold','FontSize',fontsize+1,...
    'FontAngle','italic','HorizontalAlignment','left','Rotation',0); 

% Ties-to-200s
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(3,10,:)-daily_strain_rates_combo.lon_yr_100s(3,10,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(1,:)); hold on; %QIET-SQ21
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(3,14,:)-daily_strain_rates_combo.lon_yr_100s(3,14,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(2,:)); %QIET-SQ25
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(2,10,:)-daily_strain_rates_combo.lon_yr_100s(2,10,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(3,:)); %MLOW-SQ21
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(2,11,:)-daily_strain_rates_combo.lon_yr_100s(2,11,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(4,:)); %MLOW-SQ22
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(2,14,:)-daily_strain_rates_combo.lon_yr_100s(2,14,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(5,:)); %MLOW-SQ25
% strain-rate error envelopes (any data gaps --> error estimate == 0)
tt = [strain_time; flipud(strain_time)];
ee = [squeeze(daily_strain_rates_combo.lon_yr_100s(3,10,:)+3.*daily_strain_rates_combo.delta_lon_yr_100s(3,10,:)-daily_strain_rates_combo.lon_yr_100s(3,10,ID_165));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_100s(3,10,:)-3.*daily_strain_rates_combo.delta_lon_yr_100s(3,10,:)-daily_strain_rates_combo.lon_yr_100s(3,10,ID_165)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_100s(3,14,:)+3.*daily_strain_rates_combo.delta_lon_yr_100s(3,14,:)-daily_strain_rates_combo.lon_yr_100s(3,14,ID_165));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_100s(3,14,:)-3.*daily_strain_rates_combo.delta_lon_yr_100s(3,14,:)-daily_strain_rates_combo.lon_yr_100s(3,14,ID_165)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_100s(2,10,:)+3.*daily_strain_rates_combo.delta_lon_yr_100s(2,10,:)-daily_strain_rates_combo.lon_yr_100s(2,10,ID_165));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_100s(2,10,:)-3.*daily_strain_rates_combo.delta_lon_yr_100s(2,10,:)-daily_strain_rates_combo.lon_yr_100s(2,10,ID_165)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_100s(2,11,:)+3.*daily_strain_rates_combo.delta_lon_yr_100s(2,11,:)-daily_strain_rates_combo.lon_yr_100s(2,11,ID_165));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_100s(2,11,:)-3.*daily_strain_rates_combo.delta_lon_yr_100s(2,11,:)-daily_strain_rates_combo.lon_yr_100s(2,11,ID_165)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_100s(2,14,:)+3.*daily_strain_rates_combo.delta_lon_yr_100s(2,14,:)-daily_strain_rates_combo.lon_yr_100s(2,14,ID_165));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_100s(2,14,:)-3.*daily_strain_rates_combo.delta_lon_yr_100s(2,14,:)-daily_strain_rates_combo.lon_yr_100s(2,14,ID_165)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5); hold on;
% Ties-to-200s
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(3,10,:)-daily_strain_rates_combo.lon_yr_100s(3,10,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(1,:)); hold on; %QIET-SQ21
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(3,14,:)-daily_strain_rates_combo.lon_yr_100s(3,14,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(2,:)); %QIET-SQ25
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(2,10,:)-daily_strain_rates_combo.lon_yr_100s(2,10,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(3,:)); %MLOW-SQ21
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(2,11,:)-daily_strain_rates_combo.lon_yr_100s(2,11,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(4,:)); %MLOW-SQ22
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(2,14,:)-daily_strain_rates_combo.lon_yr_100s(2,14,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(5,:)); %MLOW-SQ25

colororder('default');
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_strain(1) y_plot_strain(2)])
ylabel('$\dot{\epsilon}_{lon}(t)-\dot{\epsilon}_{lon}(t_{165})$  [ yr$^{-1}$ ]','Interpreter','latex','FontSize',fontsize)
legend(['QIET' char(8211) 'SQ21'],['QIET' char(8211) 'SQ25'],['MLOW' char(8211) 'SQ21'],['MLOW' char(8211) 'SQ22'],['MLOW' char(8211) 'SQ25'],...
    'NumColumns',3,'Location','West','EdgeColor','none','color','none');    % 'Position',[0.345 0.545 0.01 0.01],'EdgeColor','none')
set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.6,'FontSize',fontsize); 
set(gca,'ytick',y_plot_strain(1)+y_plot_strain_tick:y_plot_strain_tick:y_plot_strain(2))
grid on; box on; axe3.Layer = 'top';

%%%%%%%%%%%%%%%%%%%%%%%%% 200s (in-basin) %%%%%%%%%%%%%%%%%%%%%%%%%
axes(axe4) 
colororder({'k','k'})
yyaxis right 
for i=10
    % c_dot_err
    tt = [daily_epsilon_zz(i).t22; flipud(daily_epsilon_zz(i).t22)];
    ee = [daily_epsilon_zz(i).c_dot_delta_t_combo+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165))+(3.*daily_epsilon_zz(i).c_dot_delta_t_err_combo);...
        flipud(daily_epsilon_zz(i).c_dot_delta_t_combo+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165))-(3.*daily_epsilon_zz(i).c_dot_delta_t_err_combo))];
    nanIdx = isnan(ee); ee(nanIdx) = 0;
    patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
end
for i=11
    % c_dot_err
    tt = [daily_epsilon_zz(i).t22(1:ID_198); flipud(daily_epsilon_zz(i).t22(1:ID_198))];
    ee = [daily_epsilon_zz(i).c_dot_delta_t_combo(1:ID_198)+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165))+(3.*daily_epsilon_zz(i).c_dot_delta_t_err_combo(1:ID_198));...
        flipud(daily_epsilon_zz(i).c_dot_delta_t_combo(1:ID_198)+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165))-(3.*daily_epsilon_zz(i).c_dot_delta_t_err_combo(1:ID_198)))];
    nanIdx = isnan(ee); ee(nanIdx) = 0;
    patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
end
for i=12:16
    % c_dot_err
    tt = [daily_epsilon_zz(i).t22; flipud(daily_epsilon_zz(i).t22)];
    ee = [daily_epsilon_zz(i).c_dot_delta_t_combo+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165))+(3.*daily_epsilon_zz(i).c_dot_delta_t_err_combo);...
        flipud(daily_epsilon_zz(i).c_dot_delta_t_combo+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165))-(3.*daily_epsilon_zz(i).c_dot_delta_t_err_combo))];
    nanIdx = isnan(ee); ee(nanIdx) = 0;
    patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
end
% \dot{c}
for i=10
    plot(daily_epsilon_zz(i).t22, ...
        (daily_epsilon_zz(i).c_dot_delta_t_combo)+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165)),...
        '-','LineWidth',1.5); hold on;
end
for i=11
    plot(daily_epsilon_zz(i).t22(ID_200:end), ...
        (daily_epsilon_zz(i).c_dot_delta_t_combo(ID_200:end))+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165)),...
        '-','LineWidth',1.5); hold on;
end
for i=12:1:16
    plot(daily_epsilon_zz(i).t22, ...
        (daily_epsilon_zz(i).c_dot_delta_t_combo)+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165)),...
        '-','LineWidth',1.5); hold on;
end
colororder('default');
legend('SQ21','SQ22','SQ23','SQ24','SQ25','SQ26','SQ27','NumColumns',7,...
    'Location','NorthEast','EdgeColor','none'); 
% runoff
patch_vert = ((4*y_plot_cdot_interval)/plot_cdot_range).*(nanmean(runoff_2023_nevis(:,ID(10:16,1)),2)./10); % shrink to c_dot scale
patch_x1 = racmo_time-0.5; patch_x2 = patch_x1+1;
patch_y1 = ones(length(racmo_time),1).*y_plot_cdot(1);
patch_y2 = patch_y1 + patch_vert; 
for i=1:1:length(racmo_time)
    patch([patch_x1(i) patch_x2(i) patch_x2(i) patch_x1(i)],...
          [patch_y1(i) patch_y1(i) patch_y2(i) patch_y2(i)],...
          [0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.20); hold on;
end
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_cdot(1) y_plot_cdot(2)])
ylabel('Runoff [cm w.e.] \hspace{6mm} $\dot{c}$  [ m d$^{-1}$ ] \hspace{0mm}','Interpreter','latex','FontSize',fontsize)
set(gca,'ytick',y_plot_cdot(1)+y_plot_cdot_tick:y_plot_cdot_tick:y_plot_cdot(2),'YTickLabel',run_c_dot_label);
text(t_plot(2)-0.15,y_plot_cdot(2)-1.10,'within-basin 1150s','FontName','AvenirBold','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','right')
text(t_plot(1)-0.80, y_plot_cdot(2),'d','FontName','Helvetica','FontSize',fontsize+2);
text(t_plot(1)+0.5, y_plot_cdot(2)-0.20, 'SQ21','FontName','AvenirBold','FontSize',fontsize,'Color',[0    0.4470    0.7410],'HorizontalAlignment','center');
text(t_plot(1)+1.5, y_plot_cdot(2)-0.20, 'SQ22','FontName','AvenirBold','FontSize',fontsize,'Color',[0.8500    0.3250    0.0980],'HorizontalAlignment','center');
text(t_plot(1)+2.5, y_plot_cdot(2)-0.20, 'SQ23','FontName','AvenirBold','FontSize',fontsize,'Color',[0.9290    0.6940    0.1250],'HorizontalAlignment','center');
text(t_plot(1)+3.5, y_plot_cdot(2)-0.20, 'SQ24','FontName','AvenirBold','FontSize',fontsize,'Color',[0.4940    0.1840    0.5560],'HorizontalAlignment','center');
text(t_plot(1)+4.5, y_plot_cdot(2)-0.20, 'SQ25','FontName','AvenirBold','FontSize',fontsize,'Color',[0.4660    0.6740    0.1880],'HorizontalAlignment','center');
text(t_plot(1)+5.5, y_plot_cdot(2)-0.20, 'SQ26','FontName','AvenirBold','FontSize',fontsize,'Color',[0.3010    0.7450    0.9330],'HorizontalAlignment','center');
text(t_plot(1)+7.5, y_plot_cdot(2)-0.20, 'SQ27','FontName','AvenirBold','FontSize',fontsize,'Color',[0.6350    0.0780    0.1840],'HorizontalAlignment','center');

yyaxis left % \dot{epsilon}_{lon}
for i=3 % lakes filling and draining
    rectangle('Position',[lake_fill_dates(i) y_plot_strain(1) lake_dates(i)-lake_fill_dates(i) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',lake_blue,'EdgeColor','none');
    rectangle('Position',[lake_dates(i) y_plot_strain(1) 1 y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[1 1 0.7],'EdgeColor','none'); hold on
end
for i=1:2 % other drainages
    rectangle('Position',[lake_dates(i) y_plot_strain(1) 1 y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[0.90 0.90 0.90],'EdgeColor','none');
end
for i=4:length(lake_dates) % other drainages
    rectangle('Position',[lake_dates(i) y_plot_strain(1) 1 y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[0.90 0.90 0.90],'EdgeColor','none');
end
text(t_plot(1)+0.15,-0.035,'L2A, L2B filling...',...
         'FontName','AvenirBold','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','left')
text(lake_dates(3)+0.50, y_plot_strain(1)+0.015, {'L2A';'L2B'},'FontName','AvenirBold','FontSize',fontsize+1,...
    'FontAngle','italic','HorizontalAlignment','center','Rotation',0); 

% 200's
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_200s(10,12,:)-daily_strain_rates_combo.lon_yr_200s(10,12,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(1,:)); hold on; %SQ21-23
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_200s(11,12,:)-daily_strain_rates_combo.lon_yr_200s(11,12,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(2,:)); %SQ22-23
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_200s(14,15,:)-daily_strain_rates_combo.lon_yr_200s(14,15,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(3,:)); %SQ25-26
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_200s(14,16,:)-daily_strain_rates_combo.lon_yr_200s(14,16,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(4,:)); hold on; %SQ25-27
% strain errors
tt = [strain_time; flipud(strain_time)];
ee = [squeeze(daily_strain_rates_combo.lon_yr_200s(10,12,:)+3.*daily_strain_rates_combo.delta_lon_yr_200s(10,12,:)-daily_strain_rates_combo.lon_yr_200s(10,12,ID_165));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_200s(10,12,:)-3.*daily_strain_rates_combo.delta_lon_yr_200s(10,12,:)-daily_strain_rates_combo.lon_yr_200s(10,12,ID_165)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_200s(11,12,:)+3.*daily_strain_rates_combo.delta_lon_yr_200s(11,12,:)-daily_strain_rates_combo.lon_yr_200s(11,12,ID_165));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_200s(11,12,:)-3.*daily_strain_rates_combo.delta_lon_yr_200s(11,12,:)-daily_strain_rates_combo.lon_yr_200s(11,12,ID_165)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_200s(14,15,:)+3.*daily_strain_rates_combo.delta_lon_yr_200s(14,15,:)-daily_strain_rates_combo.lon_yr_200s(14,15,ID_165));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_200s(14,15,:)-3.*daily_strain_rates_combo.delta_lon_yr_200s(14,15,:)-daily_strain_rates_combo.lon_yr_200s(14,15,ID_165)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_200s(14,16,:)+3.*daily_strain_rates_combo.delta_lon_yr_200s(14,16,:)-daily_strain_rates_combo.lon_yr_200s(14,16,ID_165));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_200s(14,16,:)-3.*daily_strain_rates_combo.delta_lon_yr_200s(14,16,:)-daily_strain_rates_combo.lon_yr_200s(14,16,ID_165)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5); hold on;
% 200's
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_200s(10,12,:)-daily_strain_rates_combo.lon_yr_200s(10,12,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(1,:)); hold on; %SQ21-23
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_200s(11,12,:)-daily_strain_rates_combo.lon_yr_200s(11,12,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(2,:)); %SQ22-23
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_200s(14,15,:)-daily_strain_rates_combo.lon_yr_200s(14,15,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(3,:)); %SQ25-26
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_200s(14,16,:)-daily_strain_rates_combo.lon_yr_200s(14,16,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(4,:)); hold on; %SQ25-27

colororder('default');
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_strain(1) y_plot_strain(2)])
ylabel('$\dot{\epsilon}_{lon}(t)-\dot{\epsilon}_{lon}(t_{165})$  [ yr$^{-1}$ ]','Interpreter','latex','FontSize',fontsize+1)
legend(['SQ21' char(8211) '23'],['SQ22' char(8211) '23'],['SQ25' char(8211) '26'],['SQ25' char(8211) '27'],'NumColumns',2,...
'Location','West','EdgeColor','none','color','none');     % 'Position',[0.2525 0.445 0.01 0.01],'EdgeColor','none') 
set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.6,'FontSize',fontsize); 
set(gca,'ytick',y_plot_strain(1)+y_plot_strain_tick:y_plot_strain_tick:y_plot_strain(2))
grid on; box on; axe4.Layer = 'top';

%%%%%%%%%%%%%%%%%%%%%%%%% 300s (in-basin) %%%%%%%%%%%%%%%%%%%%%%%%% 
axes(axe5) 
colororder({'k','k'})
yyaxis right 
for i=17:1:22
    % c_dot_err
    tt = [daily_epsilon_zz(i).t22; flipud(daily_epsilon_zz(i).t22)];
    ee = [daily_epsilon_zz(i).c_dot_delta_t_combo+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165))+(3.*daily_epsilon_zz(i).c_dot_delta_t_err_combo);...
        flipud(daily_epsilon_zz(i).c_dot_delta_t_combo+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165))-(3.*daily_epsilon_zz(i).c_dot_delta_t_err_combo))];
    nanIdx = isnan(ee); ee(nanIdx) = 0;
    patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
end
% \dot{c}
for i=17:1:22
    % find first non-nan value for c_dot timeseries
    index_nonnan(i,1) = find((~isnan(daily_epsilon_zz(i).c_dot_delta_t_combo)),1);
    
    plot(daily_epsilon_zz(i).t22, ...
        (daily_epsilon_zz(i).c_dot_delta_t_combo)+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165)),...
        '-','LineWidth',1.5); hold on;
end
colororder('default');
legend('SQ31','SQ32','SQ34','SQ35','SQ36','SQ37','NumColumns',6); 
% runoff
patch_vert = ((4*y_plot_cdot_interval)/plot_cdot_range).*(nanmean(runoff_2023_nevis(:,ID(17:21,1)),2)./10); % shrink to c_dot scale
patch_x1 = racmo_time-0.5; patch_x2 = patch_x1+1;
patch_y1 = ones(length(racmo_time),1).*y_plot_cdot(1);
patch_y2 = patch_y1 + patch_vert; 
for i=1:1:length(racmo_time)
    patch([patch_x1(i) patch_x2(i) patch_x2(i) patch_x1(i)],...
          [patch_y1(i) patch_y1(i) patch_y2(i) patch_y2(i)],...
          [0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.20); hold on;
end
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_cdot(1) y_plot_cdot(2)])
ylabel('Runoff [cm w.e.] \hspace{6mm} $\dot{c}$  [ m d$^{-1}$ ] \hspace{0mm}','Interpreter','latex','FontSize',fontsize)
set(gca,'ytick',y_plot_cdot(1)+y_plot_cdot_tick:y_plot_cdot_tick:y_plot_cdot(2),'YTickLabel',run_c_dot_label);
text(t_plot(2)-0.15,y_plot_cdot(2)-1.10,'within-basin 1350s','FontName','AvenirBold','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','right')
text(t_plot(1)-0.80, y_plot_cdot(2),'e','FontName','Helvetica','FontSize',fontsize+2);
text(t_plot(1)+0.5, y_plot_cdot(2)-0.20, 'SQ31','FontName','AvenirBold','FontSize',fontsize,'Color',[0    0.4470    0.7410],'HorizontalAlignment','center');
text(t_plot(1)+1.5, y_plot_cdot(2)-0.20, 'SQ32','FontName','AvenirBold','FontSize',fontsize,'Color',[0.8500    0.3250    0.0980],'HorizontalAlignment','center');
text(t_plot(1)+2.5, y_plot_cdot(2)-0.20, 'SQ34','FontName','AvenirBold','FontSize',fontsize,'Color',[0.9290    0.6940    0.1250],'HorizontalAlignment','center');
text(t_plot(1)+5.5, y_plot_cdot(2)-0.20, 'SQ35','FontName','AvenirBold','FontSize',fontsize,'Color',[0.4940    0.1840    0.5560],'HorizontalAlignment','center');
text(t_plot(1)+6.5, y_plot_cdot(2)-0.20, 'SQ36','FontName','AvenirBold','FontSize',fontsize,'Color',[0.4660    0.6740    0.1880],'HorizontalAlignment','center');
text(t_plot(1)+7.5, y_plot_cdot(2)-0.20, 'SQ37','FontName','AvenirBold','FontSize',fontsize,'Color',[0.3010  0.7450  0.9330],'HorizontalAlignment','center');
      
yyaxis left % \dot{epsilon}_{lon}
for i=1:1:length(lake_dates) % other drainages
    rectangle('Position',[lake_dates(i) y_plot_strain(1) 1 y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[0.90 0.90 0.90],'EdgeColor','none');
end
for i=4 % lake filling and draining
    rectangle('Position',[lake_fill_dates(i) y_plot_strain(1) lake_dates(i)-lake_fill_dates(i) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',lake_blue,'EdgeColor','none');
    rectangle('Position',[lake_dates(i) y_plot_strain(1) 1 y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[1 1 0.7],'EdgeColor','none'); hold on
end
text(t_plot(1)+0.15,-0.035,'L3A, L3B filling...',...
         'FontName','AvenirBold','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','left')
text(lake_dates(4)+0.50, y_plot_strain(1)+0.015, {'L3A';'L3B'},'FontName','AvenirBold','FontSize',fontsize+1,...
    'FontAngle','italic','HorizontalAlignment','center','Rotation',0); 

% 300s
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_300s(17,22,:)-daily_strain_rates_combo.lon_yr_300s(17,22,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(1,:)); hold on; %SQ31-37
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_300s(18,19,:)-daily_strain_rates_combo.lon_yr_300s(18,19,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(2,:)); %SQ32-34
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_300s(19,22,:)-daily_strain_rates_combo.lon_yr_300s(19,22,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(3,:)); %SQ34-37
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_300s(20,21,:)-daily_strain_rates_combo.lon_yr_300s(20,21,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(4,:)); %SQ35-36
% strain errors
tt = [strain_time; flipud(strain_time)];
ee = [squeeze(daily_strain_rates_combo.lon_yr_300s(17,22,:)+3.*daily_strain_rates_combo.delta_lon_yr_300s(17,22,:)-daily_strain_rates_combo.lon_yr_300s(17,22,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_300s(17,22,:)-3.*daily_strain_rates_combo.delta_lon_yr_300s(17,22,:)-daily_strain_rates_combo.lon_yr_300s(17,22,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.75); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_300s(18,19,:)+3.*daily_strain_rates_combo.delta_lon_yr_300s(18,19,:)-daily_strain_rates_combo.lon_yr_300s(18,19,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_300s(18,19,:)-3.*daily_strain_rates_combo.delta_lon_yr_300s(18,19,:)-daily_strain_rates_combo.lon_yr_300s(18,19,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.75); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_300s(19,22,:)+3.*daily_strain_rates_combo.delta_lon_yr_300s(19,22,:)-daily_strain_rates_combo.lon_yr_300s(19,22,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_300s(19,22,:)-3.*daily_strain_rates_combo.delta_lon_yr_300s(19,22,:)-daily_strain_rates_combo.lon_yr_300s(19,22,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.75); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_300s(20,21,:)+3.*daily_strain_rates_combo.delta_lon_yr_300s(20,21,:)-daily_strain_rates_combo.lon_yr_300s(20,21,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_300s(20,21,:)-3.*daily_strain_rates_combo.delta_lon_yr_300s(20,21,:)-daily_strain_rates_combo.lon_yr_300s(20,21,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.75); hold on;
% 300s
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_300s(17,22,:)-daily_strain_rates_combo.lon_yr_300s(17,22,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(1,:)); hold on; %SQ31-37
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_300s(18,19,:)-daily_strain_rates_combo.lon_yr_300s(18,19,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(2,:)); %SQ32-34
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_300s(19,22,:)-daily_strain_rates_combo.lon_yr_300s(19,22,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(3,:)); %SQ34-37
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_300s(20,21,:)-daily_strain_rates_combo.lon_yr_300s(20,21,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(4,:)); %SQ35-36

colororder('default');
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_strain(1) y_plot_strain(2)])
ylabel('$\dot{\epsilon}_{lon}(t)-\dot{\epsilon}_{lon}(t_{165})$  [ yr$^{-1}$ ]','Interpreter','latex','FontSize',fontsize+1)
legend(['SQ31' char(8211) '37'],['SQ32' char(8211) '34'],['SQ34' char(8211) '37'],['SQ35' char(8211) '36'],'Location','NorthWest','NumColumns',2,...
            'Location','West','EdgeColor','none','color','none')
% legend(['SQ31' char(8211) '34'],'SQ35-36','SQ35-37','Location','NorthWest','NumColumns',3,...
%         'Location','West','EdgeColor','none');     %     'Position',[0.2075 0.090 0.01 0.01],'EdgeColor','none')
set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.6,'FontSize',fontsize); 
set(gca,'ytick',y_plot_strain(1)+y_plot_strain_tick:y_plot_strain_tick:y_plot_strain(2))
xlabel('Day of Year, 2023  [ UTC ]'); 
grid on; box on; axe5.Layer = 'top';

%% print figure 
print(gcf,'-dpng','-r300',sprintf('../../paperfigs/suppfig6.2_C6_cdot_distillations_260119.png')); 