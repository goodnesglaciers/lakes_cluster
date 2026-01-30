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
% LAS 2024-10-10: days of C2 cluster
% LAS 25-04-26: Updated FASTER locations and latest strain-rate calculations
% LAS 26-01-19: +error envelopes on strain rates and c_dot

close all; clear all

% load runoff data
load('../../nevis_lakes_cluster/racmo_station_2022_index.mat') % runoff_2022_nevis(:,ID(6,1))./10 = cm w.e. at SQ13
load('../../nevis_lakes_cluster/runoff_2022_nevis.mat');
% load station names
load('station_names.mat'); station_names = upper(station_names);

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

% Tight constraints at station-days:
tight_100s = 187:196;
tight_200s = 206:217;
tight_300s = 206:212;
tight_MLOW = 207:217; % and QIET
tight_MHIH = 207:217;

%% COMBINE TIGHT AND LOOSE RECORDS FOR CDOT
% make timeseries logical for where to use 'tight' strain rates
for i=1 % MHIH
    % first fill with loose constraints (longer sliding LS window)
    % loose daily_epsilon_zz
    daily_epsilon_zz(i).t22 = daily_epsilon_zz_loose(i).t22; % time
    daily_epsilon_zz(i).c_dot_delta_t_combo = daily_epsilon_zz_loose(i).c_dot_delta_t; % c_dot
    daily_epsilon_zz(i).c_dot_delta_t_err_combo = daily_epsilon_zz_loose(i).c_dot_delta_t_err; % c_dot_err
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

for i=2:3 % MLOW, QIET
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

for i=17:21 % 300s
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

%% plotting preliminaries 
lake_dates = [195.25, 209.75, 214.25, 209.75]; % 100s, MHIH, 200s, L3B
lake_fill_dates = [176, 184, 184, 190]; % 100s, MHIH, 200s, 300s

% SQ13 not past 200 (heightened diurnals as station tips over)
[index_200] = find((daily_epsilon_zz(6).t22 >= 200),1);
daily_epsilon_zz(6).c_dot_delta_t_combo(index_200:end) = NaN;
daily_epsilon_zz(6).u_s(index_200:end) = NaN;
daily_epsilon_zz(6).w_s(index_200:end) = NaN;
daily_epsilon_zz(6).epsilon_dot_lon(index_200:end) = NaN;
daily_epsilon_zz(6).epsilon_dot_trans(index_200:end) = NaN;
daily_epsilon_zz(6).epsilon_dot_zz(index_200:end) = NaN;
% SQ15 not past 230 (fell over)
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
% t_0 for strain_time
[ID_t0] = find(strain_time>=165.01, 1, 'first');
[ID_165] = find(daily_epsilon_zz(i).t22>=165.01, 1, 'first');

%% figure -- 2022 stack: c_dot, strain rate, runoff
figure(1); clf; 
% 5-panel stack:
set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',1.15.*[1 0 18.3 24.7]);
axe1 = axes('Position',[0.075 0.810 0.86 0.182],'Box','on','NextPlot','add','XTickLabels',[]);
axe2 = axes('Position',[0.075 0.620 0.86 0.182],'Box','on','NextPlot','add','XTickLabels',[]);
axe3 = axes('Position',[0.075 0.425 0.86 0.182],'Box','on','NextPlot','add','XTickLabels',[]);
axe4 = axes('Position',[0.075 0.230 0.86 0.182],'Box','on','NextPlot','add','XTickLabels',[]);
axe5 = axes('Position',[0.075 0.035 0.86 0.182],'Box','on','NextPlot','add');

t_plot = [206 215];

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

%%%%%%%%%%%%%%%%%%%%%%%%% 100s %%%%%%%%%%%%%%%%%%%%%%%%% 
axes(axe1) 
colororder({'k','k'})
yyaxis right % \dot{c}
% c_dot
for i=4:1:9
    % c_dot_err
    tt = [daily_epsilon_zz(i).t22; flipud(daily_epsilon_zz(i).t22)];
    ee = [daily_epsilon_zz(i).c_dot_delta_t_combo+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165))+(3.*daily_epsilon_zz(i).c_dot_delta_t_err_combo);...
        flipud(daily_epsilon_zz(i).c_dot_delta_t_combo+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165))-(3.*daily_epsilon_zz(i).c_dot_delta_t_err_combo))];
    nanIdx = isnan(ee); ee(nanIdx) = 0;
    patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
end
for i=4:1:9
    % c_dot
    plot(daily_epsilon_zz(i).t22, ...
        (daily_epsilon_zz(i).c_dot_delta_t_combo)+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165)),...
        '-','LineWidth',1.5); hold on;
end
colororder('default'); CMap = get(gca, 'ColorOrder');
% runoff
patch_vert = ((4*y_plot_cdot_interval)/plot_cdot_range).*(nanmean(runoff_2022_nevis(:,ID(4:9,1)),2)./10); % shrink to c_dot scale
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
text(t_plot(2)-0.9,y_plot_cdot(2)-0.40,'within-basin 950s','FontName','AvenirBold','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','right')
text(t_plot(1)-0.75, y_plot_cdot(2),'a','FontName','Helvetica','FontSize',fontsize+2);

text(t_plot(1)+0.5, y_plot_cdot(2)-0.20, 'SQ11','FontName','AvenirBold','FontSize',fontsize,'Color',[0    0.4470    0.7410],'HorizontalAlignment','center');
text(t_plot(1)+1.5, y_plot_cdot(2)-0.20, 'SQ12','FontName','AvenirBold','FontSize',fontsize,'Color',[0.8500    0.3250    0.0980],'HorizontalAlignment','center');
text(t_plot(1)+2.5, y_plot_cdot(2)-0.20, 'SQ13','FontName','AvenirBold','FontSize',fontsize,'Color',[0.9290    0.6940    0.1250],'HorizontalAlignment','center');
text(t_plot(1)+3.5, y_plot_cdot(2)-0.20, 'SQ14','FontName','AvenirBold','FontSize',fontsize,'Color',[0.4940    0.1840    0.5560],'HorizontalAlignment','center');
text(t_plot(1)+4.5, y_plot_cdot(2)-0.20, 'SQ15','FontName','AvenirBold','FontSize',fontsize,'Color',[0.4660    0.6740    0.1880],'HorizontalAlignment','center');
text(t_plot(1)+5.5, y_plot_cdot(2)-0.20, 'SQ16','FontName','AvenirBold','FontSize',fontsize,'Color',[0.3010    0.7450    0.9330],'HorizontalAlignment','center');

yyaxis left % \dot{epsilon}_{lon}
for i=1 % lakes filling and draining
    rectangle('Position',[lake_fill_dates(i) y_plot_strain(1) lake_dates(i)-lake_fill_dates(i) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',lake_blue,'EdgeColor','none'); 
    rectangle('Position',[lake_dates(i) y_plot_strain(1) 1 y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[1 1 0.7],'EdgeColor','none');
     text(lake_fill_dates(i)+0.25,0.35,'L1A, L1B filling...',...
         'FontName','AvenirBold','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','left')
end
for i=2:1:length(lake_dates) % other drainages
    rectangle('Position',[lake_dates(i) y_plot_strain(1) 1 y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[0.90 0.90 0.90],'EdgeColor','none');
end

% strain rates
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(4,5,:)-daily_strain_rates_combo.lon_yr_100s(4,5,ID_t0)),'d-','MarkerSize',DiamondSize,'Color',CMap(1,:)); hold on; %SQ11-12
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(5,7,:)-daily_strain_rates_combo.lon_yr_100s(5,7,ID_t0)),'d-','MarkerSize',DiamondSize,'Color',CMap(2,:)); %SQ12-14
plot(strain_time(1:1950), squeeze(daily_strain_rates_combo.lon_yr_100s(8,9,1:1950)-daily_strain_rates_combo.lon_yr_100s(8,9,ID_t0)),'d-','MarkerSize',DiamondSize,'Color',CMap(3,:)); %SQ15-16
% strain errors
patch([strain_time; flipud(strain_time)], ...
    [squeeze(daily_strain_rates_combo.lon_yr_100s(4,5,:)+3.*daily_strain_rates_combo.delta_lon_yr_100s(4,5,:)-daily_strain_rates_combo.lon_yr_100s(4,5,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_100s(4,5,:)-3.*daily_strain_rates_combo.delta_lon_yr_100s(4,5,:)-daily_strain_rates_combo.lon_yr_100s(4,5,ID_t0)))],...
    [0.5 0.5 0.5],'EdgeColor', 'none','FaceAlpha', 0.75); 
patch([strain_time; flipud(strain_time)], ...
    [squeeze(daily_strain_rates_combo.lon_yr_100s(5,7,:)+3.*daily_strain_rates_combo.delta_lon_yr_100s(5,7,:)-daily_strain_rates_combo.lon_yr_100s(5,7,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_100s(5,7,:)-3.*daily_strain_rates_combo.delta_lon_yr_100s(5,7,:)-daily_strain_rates_combo.lon_yr_100s(5,7,ID_t0)))],...
    [0.5 0.5 0.5],'EdgeColor', 'none','FaceAlpha', 0.75); 
patch([strain_time(1:1950); flipud(strain_time(1:1950))], ...
    [squeeze(daily_strain_rates_combo.lon_yr_100s(8,9,1:1950)+3.*daily_strain_rates_combo.delta_lon_yr_100s(8,9,1:1950)-daily_strain_rates_combo.lon_yr_100s(8,9,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_100s(8,9,1:1950)-3.*daily_strain_rates_combo.delta_lon_yr_100s(8,9,1:1950)-daily_strain_rates_combo.lon_yr_100s(8,9,ID_t0)))],...
    [0.5 0.5 0.5],'EdgeColor', 'none','FaceAlpha', 0.75); 
% strain rates
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(4,5,:)-daily_strain_rates_combo.lon_yr_100s(4,5,ID_t0)),'d-','MarkerSize',DiamondSize,'Color',CMap(1,:)); hold on; %SQ11-12
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_100s(5,7,:)-daily_strain_rates_combo.lon_yr_100s(5,7,ID_t0)),'d-','MarkerSize',DiamondSize,'Color',CMap(2,:)); %SQ12-14
plot(strain_time(1:2000), squeeze(daily_strain_rates_combo.lon_yr_100s(8,9,1:2000)-daily_strain_rates_combo.lon_yr_100s(8,9,ID_t0)),'d-','MarkerSize',DiamondSize,'Color',CMap(3,:)); %SQ15-16

colororder('default');
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_strain(1) y_plot_strain(2)])
ylabel('$\dot{\epsilon}_{lon}(t)-\dot{\epsilon}_{lon}(t_{165})$  [ yr$^{-1}$ ]','Interpreter','latex','FontSize',fontsize+1)
legend(['SQ11' char(8211) '12'],['SQ12' char(8211) '14'],['SQ15' char(8211) '16'],'NumColumns',3,...
    'Location','West','EdgeColor','none','color','none')
    % 'Position',[0.205 0.805 0.01 0.01],'EdgeColor','none')
set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.6,'FontSize',fontsize); 
set(gca,'ytick',y_plot_strain(1)+y_plot_strain_tick:y_plot_strain_tick:y_plot_strain(2))
grid on; box on; axe1.Layer = 'top';


%%%%%%%%%%%%%%%%%%%%%%%%% 100s-to-Ties %%%%%%%%%%%%%%%%%%%%%%%%% 
axes(axe2) 
colororder({'k','k'})
yyaxis right % \dot{c}
for i=2:3
    % c_dot_err
    tt = [daily_epsilon_zz(i).t22; flipud(daily_epsilon_zz(i).t22)];
    ee = [daily_epsilon_zz(i).c_dot_delta_t_combo+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165))+(3.*daily_epsilon_zz(i).c_dot_delta_t_err_combo);...
        flipud(daily_epsilon_zz(i).c_dot_delta_t_combo+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165))-(3.*daily_epsilon_zz(i).c_dot_delta_t_err_combo))];
    nanIdx = isnan(ee); ee(nanIdx) = 0;
    patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
end
plot(daily_epsilon_zz(3).t22, (daily_epsilon_zz(3).c_dot_delta_t_combo)+(-1.*daily_epsilon_zz(3).c_dot_delta_t_combo(ID_165)),'-','LineWidth',1.5); hold on;
plot(daily_epsilon_zz(2).t22, (daily_epsilon_zz(2).c_dot_delta_t_combo)+(-1.*daily_epsilon_zz(2).c_dot_delta_t_combo(ID_165)),'-','LineWidth',1.5); hold on;
colororder('default');
% runoff
patch_vert = ((4*y_plot_cdot_interval)/plot_cdot_range).*(nanmean(runoff_2022_nevis(:,ID(2:3,1)),2)./10); % shrink to c_dot scale
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
text(t_plot(2)-0.9,y_plot_cdot(2)-0.40,'950s-to-tiepoints','FontName','AvenirBold','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','right')
text(t_plot(1)-0.75, y_plot_cdot(2),'b','FontName','Helvetica','FontSize',fontsize+2);

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
legend(['SQ12' char(8211) 'QIET'],['SQ16' char(8211) 'QIET'],['SQ12' char(8211) 'MLOW'],'NumColumns',3,...
    'Location','West','EdgeColor','none','color','none');    % 'Position', [0.2275 0.630 0.01 0.01],'EdgeColor','none')
set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.6,'FontSize',fontsize); 
set(gca,'ytick',y_plot_strain(1)+y_plot_strain_tick:y_plot_strain_tick:y_plot_strain(2))
grid on; box on; axe2.Layer = 'top';

%%%%%%%%%%%%%%%%%%%%%%%%
axes(axe3) % Ties-to-200s
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
patch_vert = ((4*y_plot_cdot_interval)/plot_cdot_range).*(nanmean(runoff_2022_nevis(:,ID(2:3,1)),2)./10); % shrink to c_dot scale
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
text(t_plot(2)-0.9,y_plot_cdot(2)-0.40,'tiepoints-to-1150s','FontName','AvenirBold','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','right')
text(t_plot(1)-0.75, y_plot_cdot(2),'c','FontName','Helvetica','FontSize',fontsize+2);

text(t_plot(1)+0.5, y_plot_cdot(2)-0.20, 'QIET','FontName','AvenirBold','FontSize',fontsize,'Color',[0    0.4470    0.7410],'HorizontalAlignment','center');
text(t_plot(1)+1.5, y_plot_cdot(2)-0.20, 'MLOW','FontName','AvenirBold','FontSize',fontsize,'Color',[0.8500    0.3250    0.0980],'HorizontalAlignment','center');
text(t_plot(1)+2.5, y_plot_cdot(2)-0.20, 'MHIH','FontName','AvenirBold','FontSize',fontsize,'Color',[0.9290    0.6940    0.1250],'HorizontalAlignment','center');

yyaxis left % \dot{epsilon}_{lon}

for i=1:1:length(lake_dates) % other drainages
    rectangle('Position',[lake_dates(i) y_plot_strain(1) 1 y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[0.90 0.90 0.90],'EdgeColor','none');
end

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
for i=10:1:16
% c_dot_err
    tt = [daily_epsilon_zz(i).t22; flipud(daily_epsilon_zz(i).t22)];
    ee = [daily_epsilon_zz(i).c_dot_delta_t_combo+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165))+(3.*daily_epsilon_zz(i).c_dot_delta_t_err_combo);...
        flipud(daily_epsilon_zz(i).c_dot_delta_t_combo+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165))-(3.*daily_epsilon_zz(i).c_dot_delta_t_err_combo))];
    nanIdx = isnan(ee); ee(nanIdx) = 0;
    patch(tt,ee,[0.8 0.8 0.8],'EdgeColor','none','FaceAlpha',0.5); hold on;
end
% \dot{c}
for i=10:1:16
    plot(daily_epsilon_zz(i).t22, ...
        (daily_epsilon_zz(i).c_dot_delta_t_combo)+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165)),...
        '-','LineWidth',1.5); hold on;
end
colororder('default');
legend('SQ21','SQ22','SQ23','SQ24','SQ25','SQ26','SQ27','NumColumns',7,...
    'Location','NorthEast','EdgeColor','none'); 
% runoff
patch_vert = ((4*y_plot_cdot_interval)/plot_cdot_range).*(nanmean(runoff_2022_nevis(:,ID(10:16,1)),2)./10); % shrink to c_dot scale
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
text(t_plot(2)-0.90,y_plot_cdot(2)-1.1,'within-basin 1150s','FontName','AvenirBold','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','right')
text(t_plot(1)-0.75, y_plot_cdot(2),'d','FontName','Helvetica','FontSize',fontsize+2);

text(t_plot(1)+0.5, y_plot_cdot(2)-0.20, 'SQ21','FontName','AvenirBold','FontSize',fontsize,'Color',[0    0.4470    0.7410],'HorizontalAlignment','center');
text(t_plot(1)+1.5, y_plot_cdot(2)-0.20, 'SQ22','FontName','AvenirBold','FontSize',fontsize,'Color',[0.8500    0.3250    0.0980],'HorizontalAlignment','center');
text(t_plot(1)+2.5, y_plot_cdot(2)-0.20, 'SQ23','FontName','AvenirBold','FontSize',fontsize,'Color',[0.9290    0.6940    0.1250],'HorizontalAlignment','center');
text(t_plot(1)+3.5, y_plot_cdot(2)-0.20, 'SQ24','FontName','AvenirBold','FontSize',fontsize,'Color',[0.4940    0.1840    0.5560],'HorizontalAlignment','center');
text(t_plot(1)+4.5, y_plot_cdot(2)-0.20, 'SQ25','FontName','AvenirBold','FontSize',fontsize,'Color',[0.4660    0.6740    0.1880],'HorizontalAlignment','center');
text(t_plot(1)+5.5, y_plot_cdot(2)-0.20, 'SQ26','FontName','AvenirBold','FontSize',fontsize,'Color',[0.3010    0.7450    0.9330],'HorizontalAlignment','center');
text(t_plot(1)+6.5, y_plot_cdot(2)-0.20, 'SQ27','FontName','AvenirBold','FontSize',fontsize,'Color',[0.6350    0.0780    0.1840],'HorizontalAlignment','center');

yyaxis left % \dot{epsilon}_{lon}
for i=3 % lakes filling and draining
    rectangle('Position',[lake_fill_dates(i) y_plot_strain(1) lake_dates(i)-lake_fill_dates(i) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',lake_blue,'EdgeColor','none');
    rectangle('Position',[lake_dates(i) y_plot_strain(1) 1 y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[1 1 0.7],'EdgeColor','none'); hold on
     text(lake_fill_dates(i)+0.25,0.125,'L2A, L2B filling...',...
         'FontName','AvenirBold','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','left')
end
for i=1:2 % other drainages
    rectangle('Position',[lake_dates(i) y_plot_strain(1) 1 y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[0.90 0.90 0.90],'EdgeColor','none');
end
text(214+0.45,-0.025,{'L2A','L2B'},...
         'FontName','AvenirBold','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','left')

% 200's
[index_182] = find((strain_time >= 182),1);
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_200s(10,12,:)-daily_strain_rates_combo.lon_yr_200s(10,12,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(1,:)); hold on; %SQ21-23
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_200s(11,12,:)-daily_strain_rates_combo.lon_yr_200s(11,12,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(2,:)); %SQ22-23
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_200s(14,15,:)-daily_strain_rates_combo.lon_yr_200s(14,15,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(3,:)); %SQ25-26
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_200s(14,16,:)-daily_strain_rates_combo.lon_yr_200s(14,16,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(4,:)); hold on; %SQ25-27
% strain errors
tt = [strain_time(index_182:end); flipud(strain_time(index_182:end))];
ee = [squeeze(daily_strain_rates_combo.lon_yr_200s(10,12,index_182:end)+3.*daily_strain_rates_combo.delta_lon_yr_200s(10,12,index_182:end)-daily_strain_rates_combo.lon_yr_200s(10,12,ID_165));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_200s(10,12,index_182:end)-3.*daily_strain_rates_combo.delta_lon_yr_200s(10,12,index_182:end)-daily_strain_rates_combo.lon_yr_200s(10,12,ID_165)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_200s(11,12,index_182:end)+3.*daily_strain_rates_combo.delta_lon_yr_200s(11,12,index_182:end)-daily_strain_rates_combo.lon_yr_200s(11,12,ID_165));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_200s(11,12,index_182:end)-3.*daily_strain_rates_combo.delta_lon_yr_200s(11,12,index_182:end)-daily_strain_rates_combo.lon_yr_200s(11,12,ID_165)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_200s(14,15,index_182:end)+3.*daily_strain_rates_combo.delta_lon_yr_200s(14,15,index_182:end)-daily_strain_rates_combo.lon_yr_200s(14,15,ID_165));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_200s(14,15,index_182:end)-3.*daily_strain_rates_combo.delta_lon_yr_200s(14,15,index_182:end)-daily_strain_rates_combo.lon_yr_200s(14,15,ID_165)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_200s(14,16,index_182:end)+3.*daily_strain_rates_combo.delta_lon_yr_200s(14,16,index_182:end)-daily_strain_rates_combo.lon_yr_200s(14,16,ID_165));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_200s(14,16,index_182:end)-3.*daily_strain_rates_combo.delta_lon_yr_200s(14,16,index_182:end)-daily_strain_rates_combo.lon_yr_200s(14,16,ID_165)))];
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
legend(['SQ21' char(8211) '23'],['SQ22' char(8211) '23'],['SQ25' char(8211) '26'],['SQ25' char(8211) '27'],'NumColumns',4,...
'Location','SouthWest','EdgeColor','none','color','none');  % 'Position',[0.2525 0.445 0.01 0.01],'EdgeColor','none') 
set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.6,'FontSize',fontsize); 
set(gca,'ytick',y_plot_strain(1)+y_plot_strain_tick:y_plot_strain_tick:y_plot_strain(2))
grid on; box on; axe4.Layer = 'top';

%%%%%%%%%%%%%%%%%%%%%%%%% 300s (in-basin) %%%%%%%%%%%%%%%%%%%%%%%%% 
axes(axe5) 
colororder({'k','k'})
yyaxis right 
for i=17:1:21       
% c_dot_err
    tt = [daily_epsilon_zz(i).t22; flipud(daily_epsilon_zz(i).t22)];
    ee = [daily_epsilon_zz(i).c_dot_delta_t_combo+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165))+(3.*daily_epsilon_zz(i).c_dot_delta_t_err_combo);...
        flipud(daily_epsilon_zz(i).c_dot_delta_t_combo+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165))-(3.*daily_epsilon_zz(i).c_dot_delta_t_err_combo))];
    nanIdx = isnan(ee); ee(nanIdx) = 0;
    patch(tt,ee,[0.8 0.8 0.8],'EdgeColor','none','FaceAlpha',0.5); hold on;
end
% \dot{c}
for i=17:1:20
    % find first non-nan value for c_dot timeseries
    index_nonnan(i,1) = find((~isnan(daily_epsilon_zz(i).c_dot_delta_t_combo)),1);
    
    plot(daily_epsilon_zz(i).t22, ...
        (daily_epsilon_zz(i).c_dot_delta_t_combo)+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165)),...
        '-','LineWidth',1.5); hold on;
end
for i=21
    % find first non-nan value for c_dot timeseries
    index_nonnan(i,1) = find((~isnan(daily_epsilon_zz(i).c_dot_delta_t_combo)),1);
    
    plot(daily_epsilon_zz(i).t22, ...
        (daily_epsilon_zz(i).c_dot_delta_t_combo)+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(ID_165)),...
        '-','LineWidth',1.5); hold on;
end
colororder('default');
legend('SQ31','SQ34','SQ35','SQ36','SQ37','NumColumns',5); %,...
%            'Position',[0.2075 0.125 0.01 0.01])
% runoff
patch_vert = ((4*y_plot_cdot_interval)/plot_cdot_range).*(nanmean(runoff_2022_nevis(:,ID(17:21,1)),2)./10); % shrink to c_dot scale
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
text(t_plot(2)-0.9,y_plot_cdot(2)-0.40,'within-basin 1350s','FontName','AvenirBold','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','right')
text(t_plot(1)-0.75, y_plot_cdot(2),'e','FontName','Helvetica','FontSize',fontsize+2);

text(t_plot(1)+0.5, y_plot_cdot(2)-0.20, 'SQ31','FontName','AvenirBold','FontSize',fontsize,'Color',[0    0.4470    0.7410],'HorizontalAlignment','center');
text(t_plot(1)+1.5, y_plot_cdot(2)-0.20, 'SQ34','FontName','AvenirBold','FontSize',fontsize,'Color',[0.8500    0.3250    0.0980],'HorizontalAlignment','center');
text(t_plot(1)+2.5, y_plot_cdot(2)-0.20, 'SQ35','FontName','AvenirBold','FontSize',fontsize,'Color',[0.9290    0.6940    0.1250],'HorizontalAlignment','center');
text(t_plot(1)+3.5, y_plot_cdot(2)-0.20, 'SQ36','FontName','AvenirBold','FontSize',fontsize,'Color',[0.4940    0.1840    0.5560],'HorizontalAlignment','center');
text(t_plot(1)+4.5, y_plot_cdot(2)-0.20, 'SQ37','FontName','AvenirBold','FontSize',fontsize,'Color',[0.4660    0.6740    0.1880],'HorizontalAlignment','center');
        
yyaxis left % \dot{epsilon}_{lon}
for i=4 % lake filling and draining
    rectangle('Position',[lake_fill_dates(i) y_plot_strain(1) lake_dates(i)-lake_fill_dates(i) y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',lake_blue,'EdgeColor','none');
    rectangle('Position',[lake_dates(i) y_plot_strain(1) 1 y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[1 1 0.7],'EdgeColor','none'); hold on
end
for i=1:2:length(lake_dates) % other drainages
    rectangle('Position',[lake_dates(i) y_plot_strain(1) 1 y_plot_strain(2)-y_plot_strain(1)],...
        'FaceColor',[0.90 0.90 0.90],'EdgeColor','none');
end
text(210+0.05,-0.040,'L3A',...
         'FontName','AvenirBold','FontSize',fontsize+1,'FontAngle','italic','HorizontalAlignment','left')

% 300s
[index_184] = find((strain_time >= 184),1);
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_300s(17,18,:)-daily_strain_rates_combo.lon_yr_300s(17,18,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(1,:)); hold on; %SQ31-34
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_300s(19,20,:)-daily_strain_rates_combo.lon_yr_300s(19,20,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(2,:)); %SQ35-36
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_300s(19,21,:)-daily_strain_rates_combo.lon_yr_300s(19,21,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(3,:)); %SQ35-37
% strain errors
tt = [strain_time(index_184:end); flipud(strain_time(index_184:end))];
ee = [squeeze(daily_strain_rates_combo.lon_yr_300s(17,18,index_184:end)+3.*daily_strain_rates_combo.delta_lon_yr_300s(17,18,index_184:end)-daily_strain_rates_combo.lon_yr_300s(17,18,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_300s(17,18,index_184:end)-3.*daily_strain_rates_combo.delta_lon_yr_300s(17,18,index_184:end)-daily_strain_rates_combo.lon_yr_300s(17,18,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.75); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_300s(19,20,index_184:end)+3.*daily_strain_rates_combo.delta_lon_yr_300s(19,20,index_184:end)-daily_strain_rates_combo.lon_yr_300s(19,20,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_300s(19,20,index_184:end)-3.*daily_strain_rates_combo.delta_lon_yr_300s(19,20,index_184:end)-daily_strain_rates_combo.lon_yr_300s(19,20,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.75); hold on;
ee = [squeeze(daily_strain_rates_combo.lon_yr_300s(19,21,index_184:end)+3.*daily_strain_rates_combo.delta_lon_yr_300s(19,21,index_184:end)-daily_strain_rates_combo.lon_yr_300s(19,21,ID_t0));...
    flipud(squeeze(daily_strain_rates_combo.lon_yr_300s(19,21,index_184:end)-3.*daily_strain_rates_combo.delta_lon_yr_300s(19,21,index_184:end)-daily_strain_rates_combo.lon_yr_300s(19,21,ID_t0)))];
nanIdx = isnan(ee); ee(nanIdx) = 0;
patch(tt,ee,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.75); hold on;
% 300s
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_300s(17,18,:)-daily_strain_rates_combo.lon_yr_300s(17,18,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(1,:)); hold on; %SQ31-34
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_300s(19,20,:)-daily_strain_rates_combo.lon_yr_300s(19,20,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(2,:)); %SQ35-36
plot(strain_time, squeeze(daily_strain_rates_combo.lon_yr_300s(19,21,:)-daily_strain_rates_combo.lon_yr_300s(19,21,ID_165)),'d-','MarkerSize',DiamondSize,'Color',CMap(3,:)); %SQ35-37

colororder('default');
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_strain(1) y_plot_strain(2)])
ylabel('$\dot{\epsilon}_{lon}(t)-\dot{\epsilon}_{lon}(t_{165})$  [ yr$^{-1}$ ]','Interpreter','latex','FontSize',fontsize+1)
legend(['SQ31' char(8211) '34'],'SQ35-36','SQ35-37','Location','NorthWest','NumColumns',3,...
        'Location','West','EdgeColor','none','color','none');     %     'Position',[0.2075 0.090 0.01 0.01],'EdgeColor','none')
set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.6,'FontSize',fontsize); 
set(gca,'ytick',y_plot_strain(1)+y_plot_strain_tick:y_plot_strain_tick:y_plot_strain(2))
xlabel('Day of Year, 2022  [ UTC ]'); 
grid on; box on; axe5.Layer = 'top';

%% print figure
print(gcf,'-dpng','-r300',sprintf('../../paperfigs/suppfig2.2_C2_cdot_distillations_260119.png')); 

