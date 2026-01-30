%% vels ONLY plots vs. runoff at stations
% 2024-01-18 LAS: first try at a nice runoff vs. c_dot vs. drainage
% dates figure; need to understand error levels
% 2024-01-22 LAS: plotting daily averages c_dot, u_s, runoff
% 2024-03-21 LAS: error corrected in calculation of 24/48-hr strain rates
% LAS 2024-03-26: tight and loose constraints for sliding least-squares
% LAS 2025-02-23: firming things up for paper figures; tight and loose options
% LAS 2026-01-29: finalize for repository

close all; clear all

% load runoff data
load('../../nevis_lakes_cluster/racmo_station_2022_index.mat') % runoff_2022_nevis(:,ID(6,1))./10 = cm w.e. at SQ13
load('../../nevis_lakes_cluster/runoff_2022_nevis.mat');
% load station names
load station_names.mat
num_sta = length(station_names); % number of stations

%% load c_dot timeseries -- loose
load('daily_epsilon_zz_2022R_BF2_30min_BF2_UP4_sZERO_clean_w36_t12_260119.mat') % 30-min
daily_epsilon_zz_loose = daily_epsilon_zz;

% load c_dot timeseries -- tight
load('daily_epsilon_zz_2022R_BF2_30min_BF2_UP4_sZERO_clean_w18_t6_260119.mat') % 30-min
daily_epsilon_zz_tight = daily_epsilon_zz;

% % make a choice: starting with the loose
daily_epsilon_zz = daily_epsilon_zz_loose;

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
    daily_epsilon_zz(i).u_s_combo = daily_epsilon_zz_loose(i).u_s; % u_s
    daily_epsilon_zz(i).w_s_combo = daily_epsilon_zz_loose(i).w_s; % w_s
    
    for j=1:1:length(daily_epsilon_zz(i).t22)
        if daily_epsilon_zz_tight(i).t22(j,1) >= tight_MHIH(1) && daily_epsilon_zz_tight(i).t22(j,1) <= tight_MHIH(end)+1
        % assign outputs from tight constraint (shorter sliding LS window)
        daily_epsilon_zz(i).c_dot_delta_t_combo(j,1) = daily_epsilon_zz_tight(i).c_dot_delta_t(j,1); % c_dot tight
        daily_epsilon_zz(i).u_s_combo(j,1) = daily_epsilon_zz_tight(i).u_s(j,1); % u_s tight
        daily_epsilon_zz(i).w_s_combo(j,1) = daily_epsilon_zz_tight(i).w_s(j,1); % w_s tight
        
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
    daily_epsilon_zz(i).u_s_combo = daily_epsilon_zz_loose(i).u_s; % u_s
    daily_epsilon_zz(i).w_s_combo = daily_epsilon_zz_loose(i).w_s; % w_s    
    
    for j=1:1:length(daily_epsilon_zz(i).t22)
        if daily_epsilon_zz_tight(i).t22(j,1) >= tight_MLOW(1) && daily_epsilon_zz_tight(i).t22(j,1) <= tight_MLOW(end)+1
        % assign outputs from tight constraint (shorter sliding LS window)
        daily_epsilon_zz(i).c_dot_delta_t_combo(j,1) = daily_epsilon_zz_tight(i).c_dot_delta_t(j,1); % c_dot tight
        daily_epsilon_zz(i).u_s_combo(j,1) = daily_epsilon_zz_tight(i).u_s(j,1); % u_s tight
        daily_epsilon_zz(i).w_s_combo(j,1) = daily_epsilon_zz_tight(i).w_s(j,1); % w_s tight
        
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
    daily_epsilon_zz(i).u_s_combo = daily_epsilon_zz_loose(i).u_s; % u_s
    daily_epsilon_zz(i).w_s_combo = daily_epsilon_zz_loose(i).w_s; % w_s
    
    for j=1:1:length(daily_epsilon_zz(i).t22)
        if daily_epsilon_zz_tight(i).t22(j,1) >= tight_100s(1) && daily_epsilon_zz_tight(i).t22(j,1) <= tight_100s(end)+1
        % assign outputs from tight constraint
        daily_epsilon_zz(i).c_dot_delta_t_combo(j,1) = daily_epsilon_zz_tight(i).c_dot_delta_t(j,1); % c_dot tight
        daily_epsilon_zz(i).u_s_combo(j,1) = daily_epsilon_zz_tight(i).u_s(j,1); % u_s tight
        daily_epsilon_zz(i).w_s_combo(j,1) = daily_epsilon_zz_tight(i).w_s(j,1); % w_s tight
        
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
    daily_epsilon_zz(i).u_s_combo = daily_epsilon_zz_loose(i).u_s; % u_s
    daily_epsilon_zz(i).w_s_combo = daily_epsilon_zz_loose(i).w_s; % w_s
    
    for j=1:1:length(daily_epsilon_zz(i).t22)
        if daily_epsilon_zz_tight(i).t22(j,1) >= tight_200s(1) && daily_epsilon_zz_tight(i).t22(j,1) <= tight_200s(end)+1
        % assign outputs from tight constraint
        daily_epsilon_zz(i).c_dot_delta_t_combo(j,1) = daily_epsilon_zz_tight(i).c_dot_delta_t(j,1); % c_dot tight
        daily_epsilon_zz(i).u_s_combo(j,1) = daily_epsilon_zz_tight(i).u_s(j,1); % u_s tight
        daily_epsilon_zz(i).w_s_combo(j,1) = daily_epsilon_zz_tight(i).w_s(j,1); % w_s tight
        
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
    daily_epsilon_zz(i).u_s_combo = daily_epsilon_zz_loose(i).u_s; % u_s
    daily_epsilon_zz(i).w_s_combo = daily_epsilon_zz_loose(i).w_s; % w_s
    
    for j=1:1:length(daily_epsilon_zz(i).t22)
        if daily_epsilon_zz_tight(i).t22(j,1) >= tight_300s(1) && daily_epsilon_zz_tight(i).t22(j,1) <= tight_300s(end)+1
        % assign outputs from tight constraint
        daily_epsilon_zz(i).c_dot_delta_t_combo(j,1) = daily_epsilon_zz_tight(i).c_dot_delta_t(j,1); % c_dot tight
        daily_epsilon_zz(i).u_s_combo(j,1) = daily_epsilon_zz_tight(i).u_s(j,1); % u_s tight
        daily_epsilon_zz(i).w_s_combo(j,1) = daily_epsilon_zz_tight(i).w_s(j,1); % w_s tight
        
        logical_combo(i).TF(j,1) = 1; % yes, use tight constraint
        else % keep loose constraints
            logical_combo(i).TF(j,1) = 0; % no, use loose constraint 
        end
    end
end


%% plotting preliminaries 

% 2022 lake dates
lake_dates = [194.95, 209, 214.25, 210]; % 100s, MHIH, 200s, L3B
lake_fill_dates = [176, 184, 184, 190]; % 100s, MHIH, 200s, 300s

% SQ13 not past 197
[index_200] = find((daily_epsilon_zz(6).t22 >= 197),1);
daily_epsilon_zz(6).c_dot_delta_t_combo(index_200:end) = NaN;
daily_epsilon_zz(6).u_s_combo(index_200:end) = NaN;
daily_epsilon_zz(6).w_s_combo(index_200:end) = NaN;
daily_epsilon_zz(6).epsilon_dot_lon(index_200:end) = NaN;
daily_epsilon_zz(6).epsilon_dot_trans(index_200:end) = NaN;
daily_epsilon_zz(6).epsilon_dot_zz(index_200:end) = NaN;
% SQ15 not past 210
[index_210] = find((daily_epsilon_zz(8).t22 >= 210),1);
daily_epsilon_zz(8).c_dot_delta_t_combo(index_210:end) = NaN;
daily_epsilon_zz(8).u_s_combo(index_210:end) = NaN;
daily_epsilon_zz(8).w_s_combo(index_210:end) = NaN;
daily_epsilon_zz(8).epsilon_dot_lon(index_210:end) = NaN;
daily_epsilon_zz(8).epsilon_dot_trans(index_210:end) = NaN;
daily_epsilon_zz(7).epsilon_dot_trans(index_210:end) = NaN; % SQ14 uses SQ15 for yy strain rates
daily_epsilon_zz(8).epsilon_dot_zz(index_210:end) = NaN;
daily_epsilon_zz(7).epsilon_dot_zz(index_210:end) = NaN; % SQ14 uses SQ15 for yy strain rates

% t_plot = [160 230]; % time axes limits
t_plot = [175 225]; % time axes limits
% delta_t = daily_epsilon_zz(1).delta_t; % time vector
racmo_time = 1.5:1:334.5; % racmo time vector
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

% figure -- drainage dates, horizontal vels
figure(1); clf; 
set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',1.5.*[1 1 19 22]);
axe1 = axes('Position',[0.065 0.76 0.895 0.23],'Box','on','NextPlot','add','XTickLabels',[]);
axe2 = axes('Position',[0.065 0.52 0.895 0.23],'Box','on','NextPlot','add','XTickLabels',[]);
axe3 = axes('Position',[0.065 0.28 0.895 0.23],'Box','on','NextPlot','add','XTickLabels',[]);
axe4 = axes('Position',[0.065 0.04 0.895 0.23],'Box','on','NextPlot','add');

y_plot = [0 800]; % tight
y_plot_300s = [0 800];
ry_plot = [-6 6];

axes(axe1) % 100s
colororder({'k','k'})
% runoff
yyaxis right
bar(racmo_time', nanmean(runoff_2022_nevis(:,ID(4:9,1)),2)./10,...
    0.6,'FaceColor',[0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.25);
ylabel('                         Runoff at 950s  [ cm w.e. ]'); 
ylim([ry_plot(1) ry_plot(2)])
set(gca,'ytick',0:3:6);

% u_s
yyaxis left
for i=1
    rectangle('Position',[lake_fill_dates(i) y_plot(1) lake_dates(i)-lake_fill_dates(i) y_plot(2)-y_plot(1)],...
        'FaceColor',[0.97 0.97 0.97],'EdgeColor','none'); 
    rectangle('Position',[lake_dates(i) y_plot(1) 1 y_plot(2)-y_plot(1)],...
        'FaceColor',[0.90 0.90 0.90],'EdgeColor','none');
    text(lake_dates(i)+1.5,y_plot(2)-75,'L1A, L1B','FontAngle','italic')
    text(180.5,y_plot(2)-75,'950s filling...','FontAngle','italic')
end
for i=4:1:9
    plot(daily_epsilon_zz(i).t22, ...
        daily_epsilon_zz(i).u_s_combo,...
        '-','LineWidth',1.3); hold on;
end
colororder('default');
xlim([t_plot(1) t_plot(2)]); ylim([y_plot(1) y_plot(2)])
ylabel('Horizontal Velocity, $u_{s}$  [ m yr$^{-1}$ ]','Interpreter','latex','FontSize',12)
text(t_plot(1)-3.5, y_plot(2), 'a','FontName','Helvetica','FontSize',14,'FontWeight','bold')
%text(t_plot(2)-9, y_plot(2)-50,'36-hr velocities','FontName','AvenirBold','FontSize',11)
legend('SQ11','SQ12','SQ13','SQ14','SQ15','SQ16','Location','NorthWest','NumColumns',1)
set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.9,'FontSize',11,'ytick',[200:200:1000]); 
grid on; box on; axe1.Layer = 'top';

axes(axe2) % Ties
colororder({'k','k'})
% runoff
yyaxis right
bar(racmo_time', nanmean(runoff_2022_nevis(:,ID(1:3,1)),2)./10,...
    0.6,'FaceColor',[0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.25);
ylabel('                         Runoff at 1100s  [ cm w.e. ]'); 
ylim([ry_plot(1) ry_plot(2)])
set(gca,'ytick',0:3:6);

% c_dot
yyaxis left
for i=2
    rectangle('Position',[lake_fill_dates(i) y_plot(1) lake_dates(i)-lake_fill_dates(i) y_plot(2)-y_plot(1)],...
        'FaceColor',[0.97 0.97 0.97],'EdgeColor','none'); 
    rectangle('Position',[lake_dates(i) y_plot(1) 1 y_plot(2)-y_plot(1)],...
        'FaceColor',[0.90 0.90 0.90],'EdgeColor','none'); 
    text(lake_dates(i)+1.5,y_plot(2)-75,'C2.S','FontAngle','italic')
    text(lake_fill_dates(i)+0.75,y_plot(2)-75,'C2.S filling...','FontAngle','italic')
end

rectangle('Position',[lake_dates(3) y_plot(1) 1 y_plot(2)-y_plot(1)],...
        'FaceColor',goldenrod_light,'EdgeColor','none'); 
text(lake_dates(3)+2.0,y_plot(1)+350,'L2A, L2B arrivals','FontAngle','italic')

plot(daily_epsilon_zz(3).t22, daily_epsilon_zz(3).u_s_combo,'-','LineWidth',1.3); hold on;
plot(daily_epsilon_zz(2).t22, daily_epsilon_zz(2).u_s_combo,'-','LineWidth',1.3); hold on;
plot(daily_epsilon_zz(1).t22, daily_epsilon_zz(1).u_s_combo,'-','LineWidth',1.3); hold on;
colororder('default');
xlim([t_plot(1) t_plot(2)]); ylim([y_plot(1) y_plot(2)])
ylabel('Horizontal Velocity, $u_{s}$  [ m yr$^{-1}$ ]','Interpreter','latex','FontSize',12)
text(t_plot(1)-3.5, y_plot(2), 'b','FontName','Helvetica','FontSize',14,'FontWeight','bold')
legend('QIET','MLOW','MHIH','Location','NorthWest','NumColumns',1)
set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.9,'FontSize',11,'ytick',[200:200:1000]); 
grid on; box on; axe2.Layer = 'top';

axes(axe3) % 200s
colororder({'k','k'})
% runoff
yyaxis right
bar(racmo_time', nanmean(runoff_2022_nevis(:,ID(10:16,1)),2)./10,...
    0.6,'FaceColor',[0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.25);
ylabel('                         Runoff at 1150s  [ cm w.e. ]'); 
ylim([ry_plot(1) ry_plot(2)])
set(gca,'ytick',0:3:6);

% u_s
yyaxis left
for i=3
    rectangle('Position',[lake_fill_dates(i) y_plot(1) lake_dates(i)-lake_fill_dates(i) y_plot(2)-y_plot(1)],...
        'FaceColor',[0.97 0.97 0.97],'EdgeColor','none'); 
    rectangle('Position',[lake_dates(i) y_plot(1) 1 y_plot(2)-y_plot(1)],...
        'FaceColor',goldenrod_light,'EdgeColor','none');
    text(lake_dates(i)+1.5,y_plot(2)-75,'L2A, L2B','FontAngle','italic')
    text(lake_fill_dates(i)+0.5,y_plot(2)-75,'1150s filling...','FontAngle','italic')
end
for i=10:1:16
    plot(daily_epsilon_zz(i).t22, ...
        daily_epsilon_zz(i).u_s_combo,...
        '-','LineWidth',1.3); hold on;
end
colororder('default');
xlim([t_plot(1) t_plot(2)]); ylim([y_plot(1) y_plot(2)])
ylabel('Horizontal Velocity, $u_{s}$  [ m yr$^{-1}$ ]','Interpreter','latex','FontSize',12)
text(t_plot(1)-3.5, y_plot(2), 'c','FontName','Helvetica','FontSize',14,'FontWeight','bold')
legend('SQ21','SQ22','SQ23','SQ24','SQ25','SQ26','SQ27','Location','NorthWest','NumColumns',1) 
set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.9,'FontSize',11,'ytick',[200:200:1000]); 
grid on; box on; axe3.Layer = 'top';

axes(axe4) % 300s
colororder({'k','k'})
% runoff
yyaxis right
bar(racmo_time', nanmean(runoff_2022_nevis(:,ID(17:22,1)),2)./10,...
    0.6,'FaceColor',[0    0.4470    0.7410],'EdgeColor','none','FaceAlpha',0.25);
ylabel('                         Runoff at 1350s  [ cm w.e. ]'); 
ylim([ry_plot(1) ry_plot(2)])
set(gca,'ytick',0:3:6);

% u_s
yyaxis left
for i=4
    rectangle('Position',[lake_fill_dates(i) y_plot_300s(1) lake_dates(i)-lake_fill_dates(i) y_plot_300s(2)-y_plot_300s(1)],...
        'FaceColor',[0.97 0.97 0.97],'EdgeColor','none'); 
    rectangle('Position',[lake_dates(i) y_plot_300s(1) 1 y_plot_300s(2)-y_plot_300s(1)],...
        'FaceColor',plum,'EdgeColor','none'); 
    text(lake_fill_dates(i)+0.75,y_plot_300s(2)-75,'1350s filling...','FontAngle','italic')
    text(lake_dates(i)+1.5,y_plot_300s(2)-75,'L3B moulin-type drainage','FontAngle','italic')
end
for i=17:1:21
    plot(daily_epsilon_zz(i).t22, ...
        daily_epsilon_zz(i).u_s_combo,...
        '-','LineWidth',1.3); hold on;
end

colororder('default');
xlim([t_plot(1) t_plot(2)]); ylim([y_plot_300s(1) y_plot_300s(2)])
ylabel('Horizontal Velocity, $u_{s}$  [ m yr$^{-1}$ ]','Interpreter','latex','FontSize',12)
text(t_plot(1)-3.5, y_plot_300s(2), 'd','FontName','Helvetica','FontSize',14,'FontWeight','bold')
legend('SQ31','SQ34','SQ35','SQ36','SQ37','Location','NorthWest','NumColumns',1)
set(gca,'FontName','Avenir','tickdir','in','LineWidth',0.9,'FontSize',11); 
xlabel('Day of Year, 2022 [ UTC ]'); 
grid on; box on; axe4.Layer = 'top';

% %print figure
print(gcf,'-dpng','-r300',sprintf('../../paperfigs/suppfig0.1_u_s_runoff_2022_260129.png')); 
