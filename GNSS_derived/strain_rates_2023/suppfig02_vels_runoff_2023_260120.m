%% velocity vs. runoff at stations
% 2024-01-18 LAS: first try at a nice runoff vs. c_dot vs. drainage
% dates figure; need to understand error levels
% 2024-01-22 LAS: plotting daily averages c_dot, u_s, runoff
% 2024-03-21 LAS: error corrected in calculation of 24/48-hr strain rates
% LAS 2024-03-26: tight and loose constraints for sliding least-squares
% LAS 2024-07-22: update for 2023 positions
% LAS 2025-02-23: firming things up for paper figures; tight and loose
% LAS 26-01-20: error envelopes on u_s; repository

close all; clear all

%% load runoff data
load('../../nevis_lakes_cluster/racmo_station_2023_300m_index.mat') % runoff_2023_nevis(:,ID(6,1))./10 = cm w.e. at SQ13
load('../../nevis_lakes_cluster/runoff_2023_nevis_noSK_300m.mat');
runoff_2023_nevis = runoff_2023_nevis.runoff;
% load station names
load station_names_2023.mat
num_sta = length(station_names); % number of stations

%% load vel, c_dot timeseries -- loose
load('daily_strain_rates_2023R_30min_BF2_UP4_sZERO_OUT6_clean_w36_t12_260119.mat')
load('daily_epsilon_zz_2023R_30min_BF2_UP4_sZERO_OUT6_clean_w36_t12_260119.mat')
daily_epsilon_zz_loose = daily_epsilon_zz;
daily_strain_rates_loose = daily_strain_rates_2023;
strain_time = daily_strain_rates_2023.time;
station_names_short = upper(daily_strain_rates_2023.station_names); % less SQ33

% load vel, c_dot timeseries -- tight
load('daily_strain_rates_2023R_30min_BF2_UP4_sZERO_OUT6_clean_w18_t6_260119.mat')
load('daily_epsilon_zz_2023R_30min_BF2_UP4_sZERO_OUT6_clean_w18_t6_260119.mat')
daily_epsilon_zz_tight = daily_epsilon_zz;
daily_strain_rates_tight = daily_strain_rates_2023; 

% % make a choice: starting with the loose
daily_epsilon_zz = daily_epsilon_zz_loose;
daily_strain_rates = daily_strain_rates_loose;

%% COMBINE TIGHT AND LOOSE RECORDS FOR CDOT AND U_S
% make timeseries logical for where to use 'tight' strain rates

% Tight constraints at station-days 2023:
tight_100s = 186:205;
tight_MLOW = 193:205;
tight_200s = 195:205;
tight_300s = 195:205;

for i=1:3 % MHIH, MLOW, QIET
    % first fill with loose constraints (longer sliding LS window)
    % loose daily_epsilon_zz
    daily_epsilon_zz(i).t22 = daily_epsilon_zz_loose(i).t22; % time
    daily_epsilon_zz(i).c_dot_delta_t_combo = daily_epsilon_zz_loose(i).c_dot_delta_t; % c_dot
    daily_epsilon_zz(i).u_s_combo = daily_epsilon_zz_loose(i).u_s; % u_s
    daily_epsilon_zz(i).u_s_combo_sigma = squeeze(daily_strain_rates_loose.vel_lin_fit(i,:,5))'; % u_s_delta
    daily_epsilon_zz(i).w_s_combo = daily_epsilon_zz_loose(i).w_s; % w_s
    
    for j=1:1:length(daily_epsilon_zz(i).t22)
        if daily_epsilon_zz_tight(i).t22(j,1) >= tight_MLOW(1) && daily_epsilon_zz_tight(i).t22(j,1) <= tight_MLOW(end)+1
        % assign outputs from tight constraint (shorter sliding LS window)
        daily_epsilon_zz(i).c_dot_delta_t_combo(j,1) = daily_epsilon_zz_tight(i).c_dot_delta_t(j,1); % c_dot tight
        daily_epsilon_zz(i).u_s_combo(j,1) = daily_epsilon_zz_tight(i).u_s(j,1); % u_s tight
        daily_epsilon_zz(i).u_s_combo_sigma = squeeze(daily_strain_rates_tight.vel_lin_fit(i,:,5))'; % u_s_delta
        daily_epsilon_zz(i).w_s_combo(j,1) = daily_epsilon_zz_tight(i).w_s(j,1); % w_s tight
        logical_combo(i).TF(j,1) = 1; % yes, use tight constraint
        else % keep loose constraints
        logical_combo(i).TF(j,1) = 0; % no, use loose constraint 
        end
    end
end

for i=4:9 % 100s
    % first fill with loose constraints (longer sliding LS window)
    % loose daily_epsilon_zz
    daily_epsilon_zz(i).t22 = daily_epsilon_zz_loose(i).t22; % time
    daily_epsilon_zz(i).c_dot_delta_t_combo = daily_epsilon_zz_loose(i).c_dot_delta_t; % c_dot
    daily_epsilon_zz(i).u_s_combo = daily_epsilon_zz_loose(i).u_s; % u_s
    daily_epsilon_zz(i).u_s_combo_sigma = squeeze(daily_strain_rates_loose.vel_lin_fit(i,:,5))'; % u_s_delta
    daily_epsilon_zz(i).w_s_combo = daily_epsilon_zz_loose(i).w_s; % w_s
    
    for j=1:1:length(daily_epsilon_zz(i).t22)
        if daily_epsilon_zz_tight(i).t22(j,1) >= tight_100s(1) && daily_epsilon_zz_tight(i).t22(j,1) <= tight_100s(end)+1
        % assign outputs from tight constraint
        daily_epsilon_zz(i).c_dot_delta_t_combo(j,1) = daily_epsilon_zz_tight(i).c_dot_delta_t(j,1); % c_dot tight
        daily_epsilon_zz(i).u_s_combo(j,1) = daily_epsilon_zz_tight(i).u_s(j,1); % u_s tight
        daily_epsilon_zz(i).u_s_combo_sigma = squeeze(daily_strain_rates_tight.vel_lin_fit(i,:,5))'; % u_s_delta
        daily_epsilon_zz(i).w_s_combo(j,1) = daily_epsilon_zz_tight(i).w_s(j,1); % w_s tight
        logical_combo(i).TF(j,1) = 1; % yes, use tight constraint
        else % keep loose constraints
            logical_combo(i).TF(j,1) = 0; % no, use loose constraint 
        end
    end
end

for i=10:16 % 200s
    % first fill with loose constraints (longer sliding LS window)
    % loose daily_epsilon_zz
    daily_epsilon_zz(i).t22 = daily_epsilon_zz_loose(i).t22; % time
    daily_epsilon_zz(i).c_dot_delta_t_combo = daily_epsilon_zz_loose(i).c_dot_delta_t; % c_dot
    daily_epsilon_zz(i).u_s_combo = daily_epsilon_zz_loose(i).u_s; % u_s
    daily_epsilon_zz(i).u_s_combo_sigma = squeeze(daily_strain_rates_loose.vel_lin_fit(i,:,5))'; % u_s_delta
    daily_epsilon_zz(i).w_s_combo = daily_epsilon_zz_loose(i).w_s; % w_s
    
    for j=1:1:length(daily_epsilon_zz(i).t22)
        if daily_epsilon_zz_tight(i).t22(j,1) >= tight_200s(1) && daily_epsilon_zz_tight(i).t22(j,1) <= tight_200s(end)+1
        % assign outputs from tight constraint
        daily_epsilon_zz(i).c_dot_delta_t_combo(j,1) = daily_epsilon_zz_tight(i).c_dot_delta_t(j,1)'; % c_dot tight
        daily_epsilon_zz(i).u_s_combo(j,1) = daily_epsilon_zz_tight(i).u_s(j,1); % u_s tight
        daily_epsilon_zz(i).u_s_combo_sigma = squeeze(daily_strain_rates_tight.vel_lin_fit(i,:,5)); % u_s_delta
        daily_epsilon_zz(i).w_s_combo(j,1) = daily_epsilon_zz_tight(i).w_s(j,1); % w_s tight
        logical_combo(i).TF(j,1) = 1; % yes, use tight constraint
        else % keep loose constraints
            logical_combo(i).TF(j,1) = 0; % no, use loose constraint 
        end
    end
end

for i=17:num_sta-1 % 300s
    % first fill with loose constraints (longer sliding LS window)
    % loose daily_epsilon_zz
    daily_epsilon_zz(i).t22 = daily_epsilon_zz_loose(i).t22; % time
    daily_epsilon_zz(i).c_dot_delta_t_combo = daily_epsilon_zz_loose(i).c_dot_delta_t; % c_dot
    daily_epsilon_zz(i).u_s_combo = daily_epsilon_zz_loose(i).u_s; % u_s
    daily_epsilon_zz(i).u_s_combo_sigma = squeeze(daily_strain_rates_loose.vel_lin_fit(i,:,5))'; % u_s_delta
    daily_epsilon_zz(i).w_s_combo = daily_epsilon_zz_loose(i).w_s; % w_s
    
    for j=1:1:length(daily_epsilon_zz(i).t22)
        if daily_epsilon_zz_tight(i).t22(j,1) >= tight_300s(1) && daily_epsilon_zz_tight(i).t22(j,1) <= tight_300s(end)+1
        % assign outputs from tight constraint
        daily_epsilon_zz(i).c_dot_delta_t_combo(j,1) = daily_epsilon_zz_tight(i).c_dot_delta_t(j,1); % c_dot tight
        daily_epsilon_zz(i).u_s_combo(j,1) = daily_epsilon_zz_tight(i).u_s(j,1); % u_s tight
        daily_epsilon_zz(i).u_s_combo_sigma = squeeze(daily_strain_rates_tight.vel_lin_fit(i,:,5))'; % u_s_delta
        daily_epsilon_zz(i).w_s_combo(j,1) = daily_epsilon_zz_tight(i).w_s(j,1); % w_s tight
        logical_combo(i).TF(j,1) = 1; % yes, use tight constraint
        else % keep loose constraints
            logical_combo(i).TF(j,1) = 0; % no, use loose constraint 
        end
    end
end

%% plotting preliminaries 
t_plot = [175 210]; % time axes limits
delta_t = daily_epsilon_zz(1).delta_t; % delta_t velocities
racmo_time = 1.5:1:365.5; % racmo time vector
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
goldenrod_light = [1.0 1.0 0.8];
icy_plum = [203, 183, 227]./255;

% 2023 lake dates
lake_dates = [191.2, 196.6, 202.75, 200.50, 195.25]; % 100s, MHIH (IDs 150&152), 200s, 300s, IDs 154&162
lake_fill_dates = [179, 190, 183, 190, 190]; % 100s, MHIH, 200s, 300s, IDs 154&162

% % SQ35 not 170-185 in 2023
[index_170] = find((daily_epsilon_zz(20).t22 >= 170),1);
[index_185] = find((daily_epsilon_zz(20).t22 >= 185),1);
daily_epsilon_zz(20).c_dot_delta_t_combo(index_170:index_185) = NaN;
daily_epsilon_zz(20).u_s_combo(index_170:index_185) = NaN;
daily_epsilon_zz(20).w_s(index_170:index_185) = NaN;
daily_epsilon_zz(20).epsilon_dot_lon(index_170:index_185) = NaN;
daily_epsilon_zz(20).epsilon_dot_trans(index_170:index_185) = NaN;
daily_epsilon_zz(20).epsilon_dot_zz(index_170:index_185) = NaN;

%%= figure -- drainage dates, horizontal vels
figure(1); clf; 
set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',1.5.*[30 1 19 22]);
axe1 = axes('Position',[0.065 0.76 0.895 0.23],'Box','on','NextPlot','add','XTickLabels',[]);
axe2 = axes('Position',[0.065 0.52 0.895 0.23],'Box','on','NextPlot','add','XTickLabels',[]);
axe3 = axes('Position',[0.065 0.28 0.895 0.23],'Box','on','NextPlot','add','XTickLabels',[]);
axe4 = axes('Position',[0.065 0.04 0.895 0.23],'Box','on','NextPlot','add');

y_plot = [0 1000]; % tight
y_plot_300s = [0 1050];
ry_plot = [-4 6];
ry_plot_300s = [-4 6.5];

axes(axe1) % 100s
colororder({'k','k'})
% runoff
yyaxis right
bar(racmo_time', nanmean(runoff_2023_nevis(:,ID(4:9,1)),2)./10,...
    0.6,'FaceColor',[0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.20);
ylabel('                         Runoff at 950s  [ cm w.e. ]'); 
ylim([ry_plot(1) ry_plot(2)])
set(gca,'ytick',0:2:6,'FontName','Avenir','tickdir','in','LineWidth',0.9,'FontSize',11);

% u_s
yyaxis left
for i=1
    rectangle('Position',[lake_fill_dates(i) y_plot(1) lake_dates(i)-lake_fill_dates(i) y_plot(2)-y_plot(1)],...
        'FaceColor',[0.97 0.97 0.97],'EdgeColor','none'); 
    rectangle('Position',[lake_dates(i) y_plot(1) 1 y_plot(2)-y_plot(1)],...
        'FaceColor',[0.90 0.90 0.90],'EdgeColor','none');
    text(lake_dates(i)+1.25,y_plot(2)-50,'L1A-C','FontAngle','italic')
    text(180.5,y_plot(2)-50,'950s filling...','FontAngle','italic')
end 

for i=4:1:9
    % u_s
    plot(daily_epsilon_zz(i).t22, daily_epsilon_zz(i).u_s_combo,...
        '-','LineWidth',1.3); hold on;
end
% for i=4:1:9
%     % u_s_error
%     tt = [daily_epsilon_zz(i).t22; flipud(daily_epsilon_zz(i).t22)];
%     ee = [daily_epsilon_zz(i).u_s_combo+(3.*daily_epsilon_zz(i).u_s_combo_sigma);...
%         flipud(daily_epsilon_zz(i).u_s_combo-(3.*daily_epsilon_zz(i).u_s_combo_sigma))];
%     nanIdx = isnan(ee); ee(nanIdx) = 0;
%     patch(tt,ee,[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.5); hold on;
% end
for i=4:1:9
    plot(daily_epsilon_zz(i).t22, daily_epsilon_zz(i).u_s_combo,...
        '-','LineWidth',1.3); hold on;
end

colororder('default');
xlim([t_plot(1) t_plot(2)]); ylim([y_plot(1) y_plot(2)])
ylabel('Horizontal Velocity, $u_{s}$  [ m yr$^{-1}$ ]','Interpreter','latex','FontSize',12)
text(t_plot(1)-2.25, y_plot(2)+5, 'a','FontName','Helvetica','FontSize',14,'FontWeight','bold')
legend('SQ11','SQ12','SQ13','SQ14','SQ15','SQ16','Location','NorthWest','NumColumns',1)
set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.9,'FontSize',11,'ytick',200:200:1000); 
grid on; box on; axe1.Layer = 'top';

axes(axe2) % Ties
colororder({'k','k'})
% runoff
yyaxis right
bar(racmo_time', nanmean(runoff_2023_nevis(:,ID(1:3,1)),2)./10,...
    0.6,'FaceColor',[0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.20);
ylabel('                         Runoff at 1100s  [ cm w.e. ]'); 
ylim([ry_plot(1) ry_plot(2)])
set(gca,'ytick',0:2:6,'FontName','Avenir','tickdir','in','LineWidth',0.9,'FontSize',11);

% c_dot
yyaxis left
for i=2
    rectangle('Position',[lake_fill_dates(i) y_plot(1) lake_dates(i)-lake_fill_dates(i) y_plot(2)-y_plot(1)],...
        'FaceColor',[0.97 0.97 0.97],'EdgeColor','none'); 
    rectangle('Position',[lake_dates(i) y_plot(1) 1 y_plot(2)-y_plot(1)],...
        'FaceColor',[0.90 0.90 0.90],'EdgeColor','none'); 
    text(lake_dates(i)+1.25,y_plot(2)-75,{'C5.S lakes'},'FontAngle','italic')
    text(lake_fill_dates(i)+0.5,y_plot(2)-75,{'C5.S lakes';'filling...'},'FontAngle','italic')
end
rectangle('Position',[lake_dates(3) y_plot(1) 1 y_plot(2)-y_plot(1)],...
        'FaceColor',goldenrod_light,'EdgeColor','none'); 
text(lake_dates(3)+1.25,y_plot(1)+350,'L2A, L2B arrivals','FontAngle','italic')

rectangle('Position',[lake_dates(4) y_plot(1) 1 y_plot(2)-y_plot(1)],...
        'FaceColor',plum,'EdgeColor','none'); 
text(lake_dates(4)+1.5,y_plot(2)-75,'L3A, L3B arrivals','FontAngle','italic')
text(lake_dates(4)+2.0,y_plot(1)+500,'L3C, L3D arrivals','FontAngle','italic')

plot(daily_epsilon_zz(3).t22, daily_epsilon_zz(3).u_s_combo,'-','LineWidth',1.3); hold on;
plot(daily_epsilon_zz(2).t22, daily_epsilon_zz(2).u_s_combo,'-','LineWidth',1.3); hold on;
plot(daily_epsilon_zz(1).t22, daily_epsilon_zz(1).u_s_combo,'-','LineWidth',1.3); hold on;
colororder('default');
xlim([t_plot(1) t_plot(2)]); ylim([y_plot(1) y_plot(2)])
ylabel('Horizontal Velocity, $u_{s}$  [ m yr$^{-1}$ ]','Interpreter','latex','FontSize',12)
text(t_plot(1)-2.25, y_plot(2)+5, 'b','FontName','Helvetica','FontSize',14,'FontWeight','bold')
legend('QIET','MLOW','MHIH','Location','NorthWest','NumColumns',1)
set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.9,'FontSize',11,'ytick',200:200:1000); 
grid on; box on; axe2.Layer = 'top';

axes(axe3) % 200s
colororder({'k','k'})
% runoff
yyaxis right
bar(racmo_time', nanmean(runoff_2023_nevis(:,ID(10:16,1)),2)./10,...
    0.6,'FaceColor',[0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.20);
ylabel('                         Runoff at 1150s  [ cm w.e. ]'); 
ylim([ry_plot(1) ry_plot(2)])
set(gca,'ytick',0:2:6,'FontName','Avenir','tickdir','in','LineWidth',0.9,'FontSize',11);

% u_s
yyaxis left
for i=3
    rectangle('Position',[lake_fill_dates(i) y_plot(1) lake_dates(i)-lake_fill_dates(i) y_plot(2)-y_plot(1)],...
        'FaceColor',[0.97 0.97 0.97],'EdgeColor','none'); 
    rectangle('Position',[lake_dates(i) y_plot(1) 1 y_plot(2)-y_plot(1)],...
        'FaceColor',[0.90 0.90 0.90],'EdgeColor','none');
    text(lake_dates(i)+1.5,y_plot(2)-75,'L2A, L2B','FontAngle','italic')
    text(lake_fill_dates(i)+0.75,y_plot(2)-75,'1150s filling...','FontAngle','italic')
end

rectangle('Position',[lake_dates(3) y_plot(1) 1 y_plot(2)-y_plot(1)],...
        'FaceColor',goldenrod_light,'EdgeColor','none'); 
rectangle('Position',[lake_dates(4) y_plot(1) 1 y_plot(2)-y_plot(1)],...
        'FaceColor',plum,'EdgeColor','none'); 
    
for i=10:1:16
    plot(daily_epsilon_zz(i).t22, ...
        daily_epsilon_zz(i).u_s_combo,...
        '-','LineWidth',1.3); hold on;
end
colororder('default');
xlim([t_plot(1) t_plot(2)]); ylim([y_plot(1) y_plot(2)])
ylabel('Horizontal Velocity, $u_{s}$  [ m yr$^{-1}$ ]','Interpreter','latex','FontSize',12)
text(t_plot(1)-2.25, y_plot(2)+5, 'c','FontName','Helvetica','FontSize',14,'FontWeight','bold')
legend('SQ21','SQ22','SQ23','SQ24','SQ25','SQ26','SQ27','Location','NorthWest','NumColumns',1) 
set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.9,'FontSize',11,'ytick',200:200:1000); 
grid on; box on; axe3.Layer = 'top';

axes(axe4) % 300s
colororder({'k','k'})
% runoff
yyaxis right
bar(racmo_time', nanmean(runoff_2023_nevis(:,ID(17:22,1)),2)./10,...
    0.6,'FaceColor',[0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.20);
ylabel('                         Runoff at 1350s  [ cm w.e. ]'); 
ylim([ry_plot_300s(1) ry_plot_300s(2)])
set(gca,'ytick',0:2:6,'FontName','Avenir','tickdir','in','LineWidth',0.9,'FontSize',11);

% u_s
yyaxis left
for i=4
    rectangle('Position',[lake_fill_dates(i) y_plot_300s(1) lake_dates(i)-lake_fill_dates(i) y_plot_300s(2)-y_plot_300s(1)],...
        'FaceColor',[0.97 0.97 0.97],'EdgeColor','none'); 
    rectangle('Position',[lake_dates(i) y_plot_300s(1) 1 y_plot_300s(2)-y_plot_300s(1)],...
        'FaceColor',plum,'EdgeColor','none'); 
    text(lake_fill_dates(i)+0.75,y_plot_300s(2)-105,'1350s filling...','FontAngle','italic')
    text(lake_dates(i)+1.5,y_plot_300s(2)-105,'L3A, L3B','FontAngle','italic')
end
text(lake_dates(4)+2.25,y_plot(1)+350,'L3C, L3D arrivals','FontAngle','italic')
for i=17:1:22
    plot(daily_epsilon_zz(i).t22, ...
        daily_epsilon_zz(i).u_s_combo,...
        '-','LineWidth',1.3); hold on;
end
colororder('default');
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_300s(1) y_plot_300s(2)])
ylabel('Horizontal Velocity, $u_{s}$  [ m yr$^{-1}$ ]','Interpreter','latex','FontSize',12)
text(t_plot(1)-2.25, y_plot_300s(2)+5, 'd','FontName','Helvetica','FontSize',14,'FontWeight','bold')
legend('SQ31','SQ32','SQ34','SQ35','SQ36','SQ37','Location','NorthWest','NumColumns',1)
set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.9,'FontSize',11); 
xlabel('Day of Year, 2023 [ UTC ]'); 
grid on; box on; axe4.Layer = 'top';

% %print figure
print(gcf,'-dpng','-r300',sprintf('../../paperfigs/suppfig0.2_u_s_runoff_2023_175–210_260130.png')); 

% %% figure -- drainage dates, c_dot
% figure(2); clf; 
% set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',1.5.*[1 1 19 22]);
% axe1 = axes('Position',[0.06 0.76 0.895 0.23],'Box','on','NextPlot','add','XTickLabels',[]);
% axe2 = axes('Position',[0.06 0.52 0.895 0.23],'Box','on','NextPlot','add','XTickLabels',[]);
% axe3 = axes('Position',[0.06 0.28 0.895 0.23],'Box','on','NextPlot','add','XTickLabels',[]);
% axe4 = axes('Position',[0.06 0.04 0.895 0.23],'Box','on','NextPlot','add');
% 
% y_plot = [-1.5 1.5]; % tight
% y_plot = [-1.5 1.5]; % loose
% y_plot_200s = [-1.5 1.5];
% y_plot_300s = [-1.5 2.50];
% ry_plot = [-5.0286 8];
% 
% axes(axe1) % 100s
% colororder({'k','k'})
% % runoff
% yyaxis right
% bar(racmo_time', nanmean(runoff_2023_nevis(:,ID(4:9,1)),2)./10,...
%     0.6,'FaceColor',[0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.25);
% ylabel('                         Runoff at 950s  [ cm w.e. ]'); 
% ylim([ry_plot(1) ry_plot(2)])
% set(gca,'ytick',0:2:8);
% 
% % c_dot
% yyaxis left
% for i=1
%     rectangle('Position',[lake_fill_dates(i) y_plot(1) lake_dates(i)-lake_fill_dates(i) y_plot(2)-y_plot(1)],...
%         'FaceColor',[0.97 0.97 0.97],'EdgeColor','none'); 
%     rectangle('Position',[lake_dates(i) y_plot(1) 1 y_plot(2)-y_plot(1)],...
%         'FaceColor',[0.90 0.90 0.90],'EdgeColor','none');
%     text(lake_dates(i)+1.5,y_plot(2)-0.05,'L1A-C')
%     text(180.5,y_plot(2)-0.05,'10s filling')
% end
% for i=4:1:9
%     % find first non-nan value for c_dot timeseries
%     index_nonnan(i,1) = find((~isnan(daily_epsilon_zz(i).c_dot_delta_t)),1);
% 
%     plot(daily_epsilon_zz(i).t22, ...
%         (daily_epsilon_zz(i).c_dot_delta_t_combo)+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(index_nonnan(i,1))),...
%         '-','LineWidth',1.3); hold on;
% end
% colororder('default');
% xlim([t_plot(1) t_plot(2)]); ylim([y_plot(1) y_plot(2)])
% ylabel('basal uplift rate, $\dot{c}$  [ m day$^{-1}$ ]','Interpreter','latex','FontSize',12)
% %mytitle = ['$\dot{c}$  [m day$^{-1}$]']; title(mytitle,'Interpreter','latex'); 
% legend('SQ11','SQ12','SQ13','SQ14','SQ15','SQ16','Location','NorthWest','NumColumns',3)
% set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.6,'FontSize',11); 
% grid on; box on; axe1.Layer = 'top';
% 
% axes(axe2) % Ties
% colororder({'k','k'})
% % runoff
% yyaxis right
% bar(racmo_time', nanmean(runoff_2023_nevis(:,ID(1:3,1)),2)./10,...
%     0.6,'FaceColor',[0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.25);
% ylabel('               Runoff at 1100s  [ cm w.e. ]'); 
% ylim([ry_plot(1) ry_plot(2)])
% set(gca,'ytick',0:2:8);
% 
% % c_dot
% yyaxis left
% for i=2
%     rectangle('Position',[lake_fill_dates(i) y_plot(1) lake_dates(i)-lake_fill_dates(i) y_plot(2)-y_plot(1)],...
%         'FaceColor',[0.97 0.97 0.97],'EdgeColor','none'); 
%     rectangle('Position',[lake_dates(i) y_plot(1) 1 y_plot(2)-y_plot(1)],...
%         'FaceColor',[0.90 0.90 0.90],'EdgeColor','none'); 
%     text(lake_dates(i)+1.5,y_plot(2)-0.05,'150, 152')
%     text(lake_fill_dates(i)+0.75,y_plot(2)-0.05,'150, 152')
% end
% plot(daily_epsilon_zz(3).t22, (daily_epsilon_zz(3).c_dot_delta_t_combo),'-','LineWidth',1.3); hold on;
% plot(daily_epsilon_zz(2).t22, (daily_epsilon_zz(2).c_dot_delta_t_combo),'-','LineWidth',1.3); hold on;
% plot(daily_epsilon_zz(1).t22, (daily_epsilon_zz(1).c_dot_delta_t_combo),'-','LineWidth',1.3); hold on;
% 
% colororder('default');
% xlim([t_plot(1) t_plot(2)]); ylim([y_plot(1) y_plot(2)])
% ylabel('basal uplift rate, $\dot{c}$  [ m day$^{-1}$ ]','Interpreter','latex','FontSize',12)
% %mytitle = ['$\dot{c}$  [m day$^{-1}$]']; title(mytitle,'Interpreter','latex'); 
% legend('QIET','MLOW','MHIH','Location','NorthWest','NumColumns',3)
% set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.6,'FontSize',11); 
% grid on; box on; axe2.Layer = 'top';
% 
% axes(axe3) % 200s
% colororder({'k','k'})
% % runoff
% yyaxis right
% bar(racmo_time', nanmean(runoff_2023_nevis(:,ID(10:16,1)),2)./10,...
%     0.6,'FaceColor',[0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.25);
% ylabel('                         Runoff at 1150s  [ cm w.e. ]'); 
% ylim([ry_plot(1) ry_plot(2)])
% set(gca,'ytick',0:2:8);
% 
% % c_dot
% yyaxis left
% for i=3
%     rectangle('Position',[lake_fill_dates(i) y_plot_200s(1) lake_dates(i)-lake_fill_dates(i) y_plot_200s(2)-y_plot_200s(1)],...
%         'FaceColor',[0.97 0.97 0.97],'EdgeColor','none'); 
%     rectangle('Position',[lake_dates(i) y_plot_200s(1) 1 y_plot_200s(2)-y_plot_200s(1)],...
%         'FaceColor',[0.90 0.90 0.90],'EdgeColor','none');
%     text(lake_dates(i)+1.5,y_plot_200s(2)-0.05,'L2A, L2B')
%     text(lake_fill_dates(i)+0.75,y_plot_200s(2)-0.05,'20s filling')
% end
% for i=10:1:16
%     % find first non-nan value for c_dot timeseries
%     index_nonnan(i,1) = find((~isnan(daily_epsilon_zz(i).c_dot_delta_t)),1);
% 
%     plot(daily_epsilon_zz(i).t22, ...
%         (daily_epsilon_zz(i).c_dot_delta_t_combo)+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(index_nonnan(i,1))),...
%         '-','LineWidth',1.3); hold on;
% end
% colororder('default');
% xlim([t_plot(1) t_plot(2)]); ylim([y_plot_200s(1) y_plot_200s(2)])
% ylabel('basal uplift rate, $\dot{c}$  [ m day$^{-1}$ ]','Interpreter','latex','FontSize',12)
% %mytitle = ['$\dot{c}$  [m day$^{-1}$]']; title(mytitle,'Interpreter','latex'); 
% legend('SQ21','SQ22','SQ23','SQ24','SQ25','SQ26','SQ27','Location','NorthWest','NumColumns',3) 
% set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.6,'FontSize',11); 
% grid on; box on; axe3.Layer = 'top';
% 
% axes(axe4) % 300s
% colororder({'k','k'})
% % runoff
% yyaxis right
% bar(racmo_time', nanmean(runoff_2023_nevis(:,ID(17:23,1)),2)./10,...
%     0.6,'FaceColor',[0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.25);
% ylabel('                         Runoff at 1350s  [ cm w.e. ]'); 
% ylim([ry_plot(1) ry_plot(2)])
% set(gca,'ytick',0:2:8);
% 
% % c_dot
% yyaxis left
% for i=4
%     rectangle('Position',[lake_fill_dates(i) y_plot_300s(1) lake_dates(i)-lake_fill_dates(i) y_plot_300s(2)-y_plot_300s(1)],...
%         'FaceColor',[0.97 0.97 0.97],'EdgeColor','none'); 
%     rectangle('Position',[lake_dates(i) y_plot_300s(1) 1 y_plot_300s(2)-y_plot_300s(1)],...
%         'FaceColor',[0.90 0.90 0.90],'EdgeColor','none'); 
%     text(lake_fill_dates(i)+0.75,y_plot_300s(2)-0.20,'30s filling')
%     text(lake_dates(i)+1.5,y_plot_300s(2)-0.20,'L3A–D')
% end
% for i=17:1:22
%     % find first non-nan value for c_dot timeseries
%     index_nonnan(i,1) = find((~isnan(daily_epsilon_zz(i).c_dot_delta_t)),1);
% 
%     plot(daily_epsilon_zz(i).t22, ...
%         (daily_epsilon_zz(i).c_dot_delta_t_combo)+(-1.*daily_epsilon_zz(i).c_dot_delta_t_combo(index_nonnan(i,1))),...
%         '-','LineWidth',1.3); hold on;
% end
% colororder('default');
% xlim([t_plot(1) t_plot(2)]); ylim([y_plot_300s(1) y_plot_300s(2)])
% ylabel('basal uplift rate, $\dot{c}$  [ m day$^{-1}$ ]','Interpreter','latex','FontSize',12)
% %mytitle = ['$\dot{c}$  [m day$^{-1}$]']; title(mytitle,'Interpreter','latex'); 
% legend('SQ31','SQ32','SQ34','SQ35','SQ36','SQ37','Location','NorthWest','NumColumns',3)
% set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.6,'FontSize',11); 
% xlabel('Day of Year, 2023 [ UTC ]'); 
% grid on; box on; axe4.Layer = 'top';
% 
% %print figure
% % print(gcf,'-dpng','-r300',sprintf('catalogue_paperfig_2601/cdot_runoff_2023_combo_260120.png')); 
% 
% 
% %% figure -- runoff, drainage dates, vertical vels
% figure(3); clf; 
% set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',1.5.*[1 1 19 22]);
% axe1 = axes('Position',[0.06 0.76 0.895 0.23],'Box','on','NextPlot','add','XTickLabels',[]);
% axe2 = axes('Position',[0.06 0.52 0.895 0.23],'Box','on','NextPlot','add','XTickLabels',[]);
% axe3 = axes('Position',[0.06 0.28 0.895 0.23],'Box','on','NextPlot','add','XTickLabels',[]);
% axe4 = axes('Position',[0.06 0.04 0.895 0.23],'Box','on','NextPlot','add');
% 
% y_plot = [-1.00 2.50]; % tight
% y_plot_loose = [-0.50 1.50]; % loose
% y_plot_300s = [-1.25 2.0];
% 
% axes(axe1) % 100s
% colororder({'k','k'})
% % runoff
% yyaxis right
% bar(racmo_time', nanmean(runoff_2023_nevis(:,ID(4:9,1)),2)./10,...
%     0.6,'FaceColor',[0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.25);
% ylabel('Runoff at 950s  [ cm w.e. ]'); 
% ylim([ry_plot(1) ry_plot(2)]); xlim([t_plot(1) t_plot(2)]);
% set(gca,'ytick',0:2:8); 
% 
% % w_s
% yyaxis left
% for i=1
%     rectangle('Position',[lake_fill_dates(i) y_plot(1) lake_dates(i)-lake_fill_dates(i) y_plot(2)-y_plot(1)],...
%         'FaceColor',[0.97 0.97 0.97],'EdgeColor','none'); 
%     rectangle('Position',[lake_dates(i) y_plot(1) 1 y_plot(2)-y_plot(1)],...
%         'FaceColor',[0.90 0.90 0.90],'EdgeColor','none');
%     text(lake_dates(i)+1.5,y_plot(2)-0.1,'L1A-C')
%     text(180.5,y_plot(2)-0.1,'950s filling')
% end
% for i=4:1:9
%     plot(daily_epsilon_zz(i).t22, ...
%         daily_epsilon_zz(i).w_s_combo./365.25,...
%         '-','LineWidth',1.3); hold on;
% end
% colororder('default');
% xlim([t_plot(1) t_plot(2)]); ylim([y_plot(1) y_plot(2)])
% ylabel('Vertical Velocity, $w_{s}$  [ m d$^{-1}$ ]','Interpreter','latex','FontSize',12)
% text(t_plot(1)-5, y_plot(2), 'a','FontName','Helvetica','FontSize',14,'FontWeight','bold')
% legend('SQ11','SQ12','SQ13','SQ14','SQ15','SQ16','Location','NorthWest','NumColumns',3)
% set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.6,'FontSize',11); 
% grid on; box on; axe1.Layer = 'top';
% 
% axes(axe2) % Ties
% colororder({'k','k'})
% % runoff
% yyaxis right
% bar(racmo_time', nanmean(runoff_2023_nevis(:,ID(1:3,1)),2)./10,...
%     0.6,'FaceColor',[0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.25);
% ylabel('Runoff at 1100s  [ cm w.e. ]'); 
% ylim([ry_plot(1) ry_plot(2)]); xlim([t_plot(1) t_plot(2)]);
% set(gca,'ytick',0:2:8);
% 
% % c_dot
% yyaxis left
% for i=2
%     rectangle('Position',[lake_fill_dates(i) y_plot(1) lake_dates(i)-lake_fill_dates(i) y_plot(2)-y_plot(1)],...
%         'FaceColor',[0.97 0.97 0.97],'EdgeColor','none'); hold on;
%     rectangle('Position',[lake_dates(i) y_plot(1) 1 y_plot(2)-y_plot(1)],...
%         'FaceColor',[0.90 0.90 0.90],'EdgeColor','none'); 
%     text(lake_dates(i)+1.5,y_plot_loose(2)-0.06,'L150, L152')
%     text(lake_fill_dates(i)+0.75,y_plot_loose(2)-0.06,'L150, L152')
% end
% plot(daily_epsilon_zz(3).t22, daily_epsilon_zz(3).w_s_combo./365.25,'-','LineWidth',1.3); 
% plot(daily_epsilon_zz(2).t22, daily_epsilon_zz(2).w_s_combo./365.25,'-','LineWidth',1.3); 
% plot(daily_epsilon_zz(1).t22, daily_epsilon_zz(1).w_s_combo./365.25,'-','LineWidth',1.3); 
% colororder('default');
% xlim([t_plot(1) t_plot(2)]); ylim([y_plot_loose(1) y_plot_loose(2)])
% ylabel('Vertical Velocity, $w_{s}$  [ m d$^{-1}$ ]','Interpreter','latex','FontSize',12)
% text(t_plot(1)-5, y_plot_loose(2), 'b','FontName','Helvetica','FontSize',14,'FontWeight','bold')
% legend('QIET','MLOW','MHIH','Location','NorthWest','NumColumns',3)
% set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.6,'FontSize',11); 
% grid on; box on; axe2.Layer = 'top';
% 
% axes(axe3) % 200s
% colororder({'k','k'})
% % runoff
% yyaxis right
% bar(racmo_time', nanmean(runoff_2023_nevis(:,ID(10:16,1)),2)./10,...
%     0.6,'FaceColor',[0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.25);
% ylabel('Runoff at 1150s  [ cm w.e. ]'); 
% ylim([ry_plot(1) ry_plot(2)])
% set(gca,'ytick',0:2:8);
% 
% % u_s
% yyaxis left
% for i=3
%     rectangle('Position',[lake_fill_dates(i) y_plot(1) lake_dates(i)-lake_fill_dates(i) y_plot(2)-y_plot(1)],...
%         'FaceColor',[0.97 0.97 0.97],'EdgeColor','none'); 
%     rectangle('Position',[lake_dates(i) y_plot(1) 1 y_plot(2)-y_plot(1)],...
%         'FaceColor',[0.90 0.90 0.90],'EdgeColor','none');
%     text(lake_dates(i)+1.5,y_plot_loose(2)-0.1,'L2A, L2B')
%     text(lake_fill_dates(i)+0.75,y_plot_loose(2)-0.1,'20s filling')
% end
% for i=10:1:16
%     plot(daily_epsilon_zz(i).t22, ...
%         daily_epsilon_zz(i).w_s_combo./365.25,...
%         '-','LineWidth',1.3); hold on;
% end
% colororder('default');
% xlim([t_plot(1) t_plot(2)]); ylim([y_plot_loose(1) y_plot_loose(2)])
% ylabel('Vertical Velocity, $w_{s}$  [ m d$^{-1}$ ]','Interpreter','latex','FontSize',12)
% text(t_plot(1)-5, y_plot_loose(2), 'c','FontName','Helvetica','FontSize',14,'FontWeight','bold')
% legend('SQ21','SQ22','SQ23','SQ24','SQ25','SQ26','SQ27','Location','NorthWest','NumColumns',3) 
% set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.6,'FontSize',11); 
% grid on; box on; axe3.Layer = 'top';
% 
% axes(axe4) % 300s
% colororder({'k','k'})
% % runoff
% yyaxis right
% bar(racmo_time', nanmean(runoff_2023_nevis(:,ID(17:22,1)),2)./10,...
%     0.6,'FaceColor',[0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.25);
% ylabel('Runoff at 1350s  [ cm w.e. ]'); 
% ylim([ry_plot(1) ry_plot(2)])
% set(gca,'ytick',0:2:8);
% 
% % u_s
% yyaxis left
% for i=4
%     rectangle('Position',[lake_fill_dates(i) y_plot_300s(1) lake_dates(i)-lake_fill_dates(i) y_plot_300s(2)-y_plot_300s(1)],...
%         'FaceColor',[0.97 0.97 0.97],'EdgeColor','none'); 
%     rectangle('Position',[lake_dates(i) y_plot_300s(1) 1 y_plot_300s(2)-y_plot_300s(1)],...
%         'FaceColor',[0.90 0.90 0.90],'EdgeColor','none'); 
%     text(lake_fill_dates(i)+0.75,y_plot_300s(2)-0.2,'30s filling')
%     text(lake_dates(i)+1.5,y_plot_300s(2)-0.2,'L3A–D')
% end
% for i=17:1:22
%     plot(daily_epsilon_zz(i).t22, ...
%         daily_epsilon_zz(i).w_s_combo./365.25,...
%         '-','LineWidth',1.3); hold on;
% end
% colororder('default');
% xlim([t_plot(1) t_plot(2)]); ylim([y_plot_300s(1) y_plot_300s(2)])
% ylabel('Vertical Velocity, $w_{s}$  [ m d$^{-1}$ ]','Interpreter','latex','FontSize',12)
% text(t_plot(1)-5, y_plot_300s(2), 'd','FontName','Helvetica','FontSize',14,'FontWeight','bold')
% legend('SQ31','SQ32','SQ34','SQ35','SQ36','SQ37','Location','NorthWest','NumColumns',3)
% set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.6,'FontSize',11); 
% xlabel('Day of Year, 2023 [ UTC ]'); 
% grid on; box on; axe4.Layer = 'top';
% 
% % % %print figure
% % print(gcf,'-dpng','-r300',sprintf('catalogue_paperfig_2601/w_s_runoff_2023_combo_260120.png')); 
