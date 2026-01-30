%% try bed opening (c_dot) in 100s in 2023
% LAS 2024-01-11: first go; following Andrews et al. (2018) \dot\epsilon_{zz,bg}
% LAS 2024-01-15: Now with \dot\epsilon_{zz} derived from s\dot\epsilon_{xx} 
% and \dot\epsilon_{yy} of station pairs
% LAS 2024-01-16: Now with theta_b defined by 2022/150-155 \dot\epsilon_{zz} 
% (200s have uplift movement earlier than the 100s)
% LAS 2024-01-17: Now with all stations, except SQ33 (slide-away). 24-hr
% and 48-hr windows for linear regressions to calculate vertical
% velocities, need at least 24 hours of points to keep c_dot estimate.
% LAS 2024-03-21: FOUND AN ERROR -- had been using horizontal velocities
% for vertical velocities in the calculation of strain rates by switching 
% the order from NEU to ENU !!!
% LAS 2024-03-26: tight and loose constraints for sliding least-squares
% LAS 2024-07-10: adjust for 2023 season
% LAS 2024-08-16: double check for sign errors in 2023 data, testing
% different bias-flag thresholds
% LAS 2024-10-30: attempt using ATM and RMS stats 
% LAS 2025-04: remove gnarly edges, 30-min velocities
% LAS 2026-01-13: keep track of data-quality rejected proportions

%% variables needed for each GNSS station
% w_s = vertical speed [m/yr] 
% w_s_bg = vertical speed background (DOY 150–160) [m/yr]
% u_s = horizontal velocity at surface [m/yr]
% u_s_bg = horizontal velocity at surface background (DOY 150-160) [m/yr]
% (assume equal to u_b_bg)
% u_b = horizontal velocity at bed [m/yr]
% u_b_bg = horizontal velocity at bed background (DOY 150-160) [m/yr]
% z_s = GNSS vertical displacement [m]
% theta_b = local bed slope [deg]
% H = ice thickness [m]
% epsilon_dot_xx = along-flow strain rate [1/yr]
% epsilon_dot_yy = across-flow strain rate [1/yr]
% epsilon_dot_zz = vertical strain rate [1/yr]
% delta_t = time step [yr]

% solve for
% c_dot = basal opening rate [m/yr]
% c_dot_delta_t = basal opening [m] (integrates from start of timeseries)

clear all; close all
%% load 2022 meltseason positions
load A2022_unfixedBF2_atm_datarejected_260112.mat % only data with <2 unfixed bias flags
station_names = {A2022.name}; sites = station_names;
num_sta = length(station_names)-1; % number of stations; omit SQ33
station_names_flip = vertcat(station_names(22),station_names(21),station_names(20),station_names(19),...
    station_names(18),station_names(17),station_names(16),...
    station_names(15),station_names(14),station_names(13),...
    station_names(12),station_names(11),station_names(10),...
    station_names(9),station_names(8),station_names(7),...
    station_names(6),station_names(5),station_names(4),...
    station_names(3),station_names(2),station_names(1));
save station_names.mat station_names
station_names = horzcat(station_names(1:17), station_names(20:end)); % omit SQ33

%% load north lake geographic files and morlighem bed relative to North Lake
% XY relative to moulin M1 in Stevens et al. (2015)
origin = [68.72, -49.53]; % M1 moulin
% polarstereo conversion needed values
radius=6378137.0;    eccen=0.08181919;    lat_true=70;    lon_posy=-45;
[moulin_x,moulin_y] = polarstereo_fwd(origin(1),origin(2),radius,eccen,lat_true,lon_posy); % [ m ]
moulin_x_km = moulin_x./1e3; moulin_y_km = moulin_y./1e3; % UTM in KM of moulin
% load station names
station_names_2022 = load('station_names.mat');
num_sta_22 = length(station_names_2022.station_names)-2; % number of stations
station_names_2022_short = upper(vertcat(station_names_2022.station_names(1,1:17)',...
    station_names_2022.station_names(1,20:22)')); % omit SQ32, SQ33 (#18, #19)
% load site locations 2022
load('apcoords_lle_2022.mat');
    lats=apcoords_lle_2022(1,:); lons=apcoords_lle_2022(2,:)-360; hs=apcoords_lle_2022(3,:); llh=[lats; lons; hs];
    [xy_sta_22(:,1), xy_sta_22(:,2)] = polarstereo_fwd(lats,lons,radius,eccen,lat_true,lon_posy); % [ m ]
    xy_sta_22(:,1) = (xy_sta_22(:,1)-moulin_x)./1e3; % [ km ] 
    xy_sta_22(:,2) = (xy_sta_22(:,2)-moulin_y)./1e3; % [ km ] 
    xy_sta_22_short = vertcat(xy_sta_22(1:17,:), xy_sta_22(19:end,:)); % omit SQ33 (#18)
save polarstereo_stations_2022_short.mat xy_sta_22_short

%% load background variables (H, bed slope, and background speeds)
% w_s_bg = vertical speed background (DOY 140–160) [m/yr]
% u_s_bg = horizontal velocity at surface background (DOY 150–160) [m/yr]
% (assume equal to u_b_bg)
% u_b_bg = horizontal velocity at bed background (DOY 150–160) [m/yr]
% theta_b = local bed slope [deg]
% H = ice thickness [m]
load('environs_2022_150_160_260119.mat') % 2022/150-160 = background date window
w_s_bg = environs.w_s_bg_flow_vel_myr; % [ m/yr]
u_s_bg = environs.u_s_bg_flow_vel_myr; % [ m/yr]
u_b_bg = u_s_bg; % [ m/yr ] assume u_b_bg = u_s_bg
theta_b = environs.bed_slope_GNSS_deg; % [ degrees ] (positive = downhill with flow)
theta_b_err = environs.bed_slope_GNSS_deg_err; % [ degrees ] (all positive)
H = environs.thickness_m; % [ m ]
H_err = environs.thickness_error_m; % [ m ]

% omit SQ33 (#18)
w_s_bg = vertcat(w_s_bg(1:17), w_s_bg(19:end)); % omit SQ33 (#18)
u_s_bg = vertcat(u_s_bg(1:17), u_s_bg(19:end)); % omit SQ33
u_b_bg = vertcat(u_b_bg(1:17), u_b_bg(19:end)); % omit SQ33
theta_b = vertcat(theta_b(1:17), theta_b(19:end)); % omit SQ33
theta_b_err = vertcat(theta_b_err(1:17), theta_b_err(19:end)); % omit SQ33
max(abs(theta_b)) % small-angle assumption OK bc theta_b<10 degrees (for error propagation)
H = vertcat(H(1:17), H(19:end)); % omit SQ33
H_err = vertcat(H_err(1:17), H_err(19:end)); % omit SQ33

%% solve for background epsilon_dot_zz_bg [ 1/yr ]
% assume u_b_bg = u_s_bg
%%%% CHECK SIGN CONVENTION FOR BED SLOPE! %%%%
epsilon_dot_zz_bg = (w_s_bg - (u_b_bg.*tand(-1.*theta_b)))./H; % [ 1/yr ]
% mean(abs(epsilon_dot_zz_bg)) = 0.004 yr^-1

% figure(1); clf
% set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',1.2.*[5 2 22 12]);
% axe1 = axes('Position',[0.08 0.1 0.9 0.85],'Box','on','NextPlot','add','XTickLabels',[]);
% axes(axe1)
% plot(epsilon_dot_zz_bg,'ks'); hold on;
% set(gca,'xtick',1:1:21,'xticklabel',station_names);
% %for i=1:1:22; text(i+0.2,epsilon_dot_zz_bg(i,1),sprintf('%5.4f',epsilon_dot_zz_bg(i,1))); end
% grid on; xlim([0 22]); ylim([-0.01 0.01])
% title('thickness-integrated vertical strain rate during 2022/150-155  [ yr^{-1} ]')
% ylabel('$\dot{\epsilon}_{zz,bg}$  [ yr$^{-1}$ ]','Interpreter','latex','FontSize',14)

% print figure
% cd catalogue_c_dot/100s/
% print(gcf,'-dpng','-r500',sprintf('epsilonzz_2022_150_170_20240114.png')); 
% cd ../../

%% time preliminaries
time_start = 150; time_end = 230;
interptime_15sec = 1.736111100001381e-04;
ti22_15sec = (time_start:interptime_15sec:time_end)'; % 15-second time vector

%% remove SQ32, SQ33 from the analysis 
A2023 = A2022; % switch to 2023-code variable names
for i=1:1:17
    % time indices
    idx_time(i,1) = find((time_start-A2023(i).neu(:,2))<=0,1); % start
    idx_time(i,2) = find((time_end+2-A2023(i).neu(:,2))<=0,1); % stop

    % timeseries
    B2023s(i).neu(:,1) = A2023(i).neu(idx_time(i,1):idx_time(i,2),2); 
    B2023s(i).neu(:,2) = A2023(i).neu(idx_time(i,1):idx_time(i,2),3);
    B2023s(i).neu(:,3) = A2023(i).neu(idx_time(i,1):idx_time(i,2),5);
    B2023s(i).neu(:,4) = A2023(i).neu(idx_time(i,1):idx_time(i,2),7);
    B2023s(i).neu(:,5) = A2023(i).neu(idx_time(i,1):idx_time(i,2),4); % 1-sigma Track error north
    B2023s(i).neu(:,6) = A2023(i).neu(idx_time(i,1):idx_time(i,2),6); % 1-sigma Track error east
    B2023s(i).neu(:,7) = A2023(i).neu(idx_time(i,1):idx_time(i,2),8); % 1-sigma Track error up

    % rejected data
    A2023(i).data_rejected(A2023(i).data_rejected(:,1) == 2023, :) = [];
    % time indices
    rdx_time(i,1) = find((A2023(i).data_rejected(:,1)-time_start)>=0,1); % start
    rdx_time(i,2) = find((A2023(i).data_rejected(:,1)-time_end-2)<=0,1,'last'); % stop
    % rejected data due to bias flags > 2
    R2022(i).data_rejected_BF = A2023(i).data_rejected(rdx_time(i,1):rdx_time(i,2),:); % DOY, data points total, points rejected, points usable, proportion usable
end

for i=20:1:23 % omit SQ32, SQ33 (rows 18, 19)
    % time indices
    idx_time(i-2,1) = find((time_start-A2023(i).neu(:,2))<=0,1); % start
    idx_time(i-2,2) = find((time_end+2-A2023(i).neu(:,2))<=0,1); % stop

    % timeseries
    B2023s(i-2).neu(:,1) = A2023(i).neu(idx_time(i-2,1):idx_time(i-2,2),2); 
    B2023s(i-2).neu(:,2) = A2023(i).neu(idx_time(i-2,1):idx_time(i-2,2),3);
    B2023s(i-2).neu(:,3) = A2023(i).neu(idx_time(i-2,1):idx_time(i-2,2),5);
    B2023s(i-2).neu(:,4) = A2023(i).neu(idx_time(i-2,1):idx_time(i-2,2),7);
    B2023s(i-2).neu(:,5) = A2023(i).neu(idx_time(i-2,1):idx_time(i-2,2),4); % 1-sigma Track error north
    B2023s(i-2).neu(:,6) = A2023(i).neu(idx_time(i-2,1):idx_time(i-2,2),6); % 1-sigma Track error east
    B2023s(i-2).neu(:,7) = A2023(i).neu(idx_time(i-2,1):idx_time(i-2,2),8); % 1-sigma Track error up

    % rejected data
    A2023(i).data_rejected(A2023(i).data_rejected(:,1) == 2023, :) = [];
    % time indices
    rdx_time(i-2,1) = find((A2023(i).data_rejected(:,1)-time_start)>=0,1); % start
    rdx_time(i-2,2) = find((A2023(i).data_rejected(:,1)-time_end-2)<=0,1,'last'); % stop
    % rejected data due to bias flags > 2
    R2022(i-2).data_rejected_BF = A2023(i).data_rejected(rdx_time(i-2,1):rdx_time(i-2,2),:); % DOY, data points total, points rejected, points usable, proportion usable
end

%figure; plot(B2023s(3).neu(:,1), B2023s(3).neu(:,3),'.');
%figure; plot(B2023s(18).neu(:,1), B2023s(18).neu(:,4),'.');

%% remove outliers more than 3 sigma from mean
% identify and remove outliers over X hours
span_outliers = (4*60)*6;  % window width [points] [6 hours]
% 12 hours or longer cuts out the lake drainage for some stations !!
days_EXEMPT = [195, 195, 209, 210, 214, 215]; % DOY of lake drainages

for i=1 % MHIH
    % station names
    B2023(i).site = station_names(i);
    % outlier identification 
    [U2023(i).n22, TFrm(i).n22] = rmoutliers(B2023s(i).neu(:,2),'movmedian',span_outliers); % N
    [U2023(i).e22, TFrm(i).e22] = rmoutliers(B2023s(i).neu(:,3),'movmedian',span_outliers); % E
    [U2023(i).u22, TFrm(i).u22] = rmoutliers(B2023s(i).neu(:,4),'movmedian',span_outliers); % U
    % keep non-outliers: any direction of outlier makes the time point an outlier
    TFrm(i).enu_initial = TFrm(i).n22 + TFrm(i).e22 + TFrm(i).u22; 
    % don't worry about lake-drainage days
    NO_DRAIN = round(B2023s(i).neu(:,1)) ~= days_EXEMPT(3:4);
    NO_DRAIN_ANY = NO_DRAIN(:,1).*NO_DRAIN(:,2); % lake drainages are 0's
    TFrm(i).enu = TFrm(i).enu_initial.*NO_DRAIN_ANY;
    % remove outliers
    B2023(i).neu_iu(:,1) = B2023s(i).neu(~TFrm(i).enu,1); % time
    B2023(i).neu_iu(:,2) = B2023s(i).neu(~TFrm(i).enu,2); % north
    B2023(i).neu_iu(:,3) = B2023s(i).neu(~TFrm(i).enu,3); % east
    B2023(i).neu_iu(:,4) = B2023s(i).neu(~TFrm(i).enu,4); % up
    B2023(i).neu_iu(:,5) = B2023s(i).neu(~TFrm(i).enu,5); % north 1-sigma
    B2023(i).neu_iu(:,6) = B2023s(i).neu(~TFrm(i).enu,6); % east 1-sigma
    B2023(i).neu_iu(:,7) = B2023s(i).neu(~TFrm(i).enu,7); % up 1-sigma
    % rejected data due to outliers
    R2022(i).data_rejected_outliers(1,1:4) = [length(B2023s(i).neu(:,2)), sum(TFrm(i).enu), ...
        length(B2023s(i).neu(:,2))-sum(TFrm(i).enu), sum(TFrm(i).enu)./length(B2023s(i).neu(:,2))]; % data points total, points rejected, points usable, proportion unusable 
end

for i=2:3 % MLOW, QIET
    % station names
    B2023(i).site = station_names(i);
    % outlier identification 
    [U2023(i).n22, TFrm(i).n22] = rmoutliers(B2023s(i).neu(:,2),'movmedian',span_outliers); % N
    [U2023(i).e22, TFrm(i).e22] = rmoutliers(B2023s(i).neu(:,3),'movmedian',span_outliers); % E
    [U2023(i).u22, TFrm(i).u22] = rmoutliers(B2023s(i).neu(:,4),'movmedian',span_outliers); % U
    % keep non-outliers: any direction of outlier makes the time point an outlier
    TFrm(i).enu = TFrm(i).n22 + TFrm(i).e22 + TFrm(i).u22; 
    % remove outliers
    B2023(i).neu_iu(:,1) = B2023s(i).neu(~TFrm(i).enu,1); % time
    B2023(i).neu_iu(:,2) = B2023s(i).neu(~TFrm(i).enu,2); % north
    B2023(i).neu_iu(:,3) = B2023s(i).neu(~TFrm(i).enu,3); % east
    B2023(i).neu_iu(:,4) = B2023s(i).neu(~TFrm(i).enu,4); % up
    B2023(i).neu_iu(:,5) = B2023s(i).neu(~TFrm(i).enu,5); % north 1-sigma
    B2023(i).neu_iu(:,6) = B2023s(i).neu(~TFrm(i).enu,6); % east 1-sigma
    B2023(i).neu_iu(:,7) = B2023s(i).neu(~TFrm(i).enu,7); % up 1-sigma
    % rejected data due to outliers
    R2022(i).data_rejected_outliers(1,1:4) = [length(B2023s(i).neu(:,2)), sum(TFrm(i).enu), ...
        length(B2023s(i).neu(:,2))-sum(TFrm(i).enu), sum(TFrm(i).enu)./length(B2023s(i).neu(:,2))]; % data points total, points rejected, points usable, proportion unusable 
end

for i=4:9 % 10s
    % station names
    B2023(i).site = station_names(i);
    % outlier identification 
    [U2023(i).n22, TFrm(i).n22] = rmoutliers(B2023s(i).neu(:,2),'movmedian',span_outliers); % N
    [U2023(i).e22, TFrm(i).e22] = rmoutliers(B2023s(i).neu(:,3),'movmedian',span_outliers); % E
    [U2023(i).u22, TFrm(i).u22] = rmoutliers(B2023s(i).neu(:,4),'movmedian',span_outliers); % U
    % keep non-outliers: any direction of outlier makes the time point an outlier
    TFrm(i).enu_initial = TFrm(i).n22 + TFrm(i).e22 + TFrm(i).u22; 
    % don't worry about lake-drainage days
    NO_DRAIN = round(B2023s(i).neu(:,1)) ~= days_EXEMPT(1);
    NO_DRAIN_ANY = NO_DRAIN(:,1); % lake drainages are 0's
    TFrm(i).enu = TFrm(i).enu_initial.*NO_DRAIN_ANY;
    % remove outliers
    B2023(i).neu_iu(:,1) = B2023s(i).neu(~TFrm(i).enu,1); % time
    B2023(i).neu_iu(:,2) = B2023s(i).neu(~TFrm(i).enu,2); % north
    B2023(i).neu_iu(:,3) = B2023s(i).neu(~TFrm(i).enu,3); % east
    B2023(i).neu_iu(:,4) = B2023s(i).neu(~TFrm(i).enu,4); % up
    B2023(i).neu_iu(:,5) = B2023s(i).neu(~TFrm(i).enu,5); % north 1-sigma
    B2023(i).neu_iu(:,6) = B2023s(i).neu(~TFrm(i).enu,6); % east 1-sigma
    B2023(i).neu_iu(:,7) = B2023s(i).neu(~TFrm(i).enu,7); % up 1-sigma
    % rejected data due to outliers
    R2022(i).data_rejected_outliers(1,1:4) = [length(B2023s(i).neu(:,2)), sum(TFrm(i).enu), ...
        length(B2023s(i).neu(:,2))-sum(TFrm(i).enu), sum(TFrm(i).enu)./length(B2023s(i).neu(:,2))]; % data points total, points rejected, points usable, proportion unusable 
end

for i=10:16 % 20s
    % station names
    B2023(i).site = station_names(i);
    % outlier identification 
    [U2023(i).n22, TFrm(i).n22] = rmoutliers(B2023s(i).neu(:,2),'movmedian',span_outliers); % N
    [U2023(i).e22, TFrm(i).e22] = rmoutliers(B2023s(i).neu(:,3),'movmedian',span_outliers); % E
    [U2023(i).u22, TFrm(i).u22] = rmoutliers(B2023s(i).neu(:,4),'movmedian',span_outliers); % U
    % keep non-outliers: any direction of outlier makes the time point an outlier
    TFrm(i).enu_initial = TFrm(i).n22 + TFrm(i).e22 + TFrm(i).u22; 
    % don't worry about lake-drainage days
    NO_DRAIN = round(B2023s(i).neu(:,1)) ~= days_EXEMPT(5:6);
    NO_DRAIN_ANY = NO_DRAIN(:,1).*NO_DRAIN(:,2); % lake drainages are 0's
    TFrm(i).enu = TFrm(i).enu_initial.*NO_DRAIN_ANY;
    % remove outliers
    B2023(i).neu_iu(:,1) = B2023s(i).neu(~TFrm(i).enu,1); % time
    B2023(i).neu_iu(:,2) = B2023s(i).neu(~TFrm(i).enu,2); % north
    B2023(i).neu_iu(:,3) = B2023s(i).neu(~TFrm(i).enu,3); % east
    B2023(i).neu_iu(:,4) = B2023s(i).neu(~TFrm(i).enu,4); % up
    B2023(i).neu_iu(:,5) = B2023s(i).neu(~TFrm(i).enu,5); % north 1-sigma
    B2023(i).neu_iu(:,6) = B2023s(i).neu(~TFrm(i).enu,6); % east 1-sigma
    B2023(i).neu_iu(:,7) = B2023s(i).neu(~TFrm(i).enu,7); % up 1-sigma
    % rejected data due to outliers
    R2022(i).data_rejected_outliers(1,1:4) = [length(B2023s(i).neu(:,2)), sum(TFrm(i).enu), ...
        length(B2023s(i).neu(:,2))-sum(TFrm(i).enu), sum(TFrm(i).enu)./length(B2023s(i).neu(:,2))]; % data points total, points rejected, points usable, proportion unusable 
end

for i=17:21 %  SQ31, SQ34-37
    % station names
    B2023(i).site = station_names(i);
    % outlier identification 
    [U2023(i).n22, TFrm(i).n22] = rmoutliers(B2023s(i).neu(:,2),'movmedian',span_outliers); % N
    [U2023(i).e22, TFrm(i).e22] = rmoutliers(B2023s(i).neu(:,3),'movmedian',span_outliers); % E
    [U2023(i).u22, TFrm(i).u22] = rmoutliers(B2023s(i).neu(:,4),'movmedian',span_outliers); % U
    % keep non-outliers: any direction of outlier makes the time point an outlier
    TFrm(i).enu_initial = TFrm(i).n22 + TFrm(i).e22 + TFrm(i).u22; 
    % don't worry about lake-drainage days
    NO_DRAIN = round(B2023s(i).neu(:,1)) ~= days_EXEMPT(3:4);
    NO_DRAIN_ANY = NO_DRAIN(:,1).*NO_DRAIN(:,2); % lake drainages are 0's
    TFrm(i).enu = TFrm(i).enu_initial.*NO_DRAIN_ANY;
    % remove outliers
    B2023(i).neu_iu(:,1) = B2023s(i).neu(~TFrm(i).enu,1); % time
    B2023(i).neu_iu(:,2) = B2023s(i).neu(~TFrm(i).enu,2); % north
    B2023(i).neu_iu(:,3) = B2023s(i).neu(~TFrm(i).enu,3); % east
    B2023(i).neu_iu(:,4) = B2023s(i).neu(~TFrm(i).enu,4); % up
    B2023(i).neu_iu(:,5) = B2023s(i).neu(~TFrm(i).enu,5); % north 1-sigma
    B2023(i).neu_iu(:,6) = B2023s(i).neu(~TFrm(i).enu,6); % east 1-sigma
    B2023(i).neu_iu(:,7) = B2023s(i).neu(~TFrm(i).enu,7); % up 1-sigma
    % rejected data due to outliers
    R2022(i).data_rejected_outliers(1,1:4) = [length(B2023s(i).neu(:,2)), sum(TFrm(i).enu), ...
        length(B2023s(i).neu(:,2))-sum(TFrm(i).enu), sum(TFrm(i).enu)./length(B2023s(i).neu(:,2))]; % data points total, points rejected, points usable, proportion unusable 
end

num_sta = length([B2023.site])

%figure; plot(B2023(18).neu_iu(:,1), B2023(18).neu_iu(:,4),'.');
%figure; plot(B2023(1).neu_iu(:,1), B2023(1).neu_iu(:,4),'.');
%figure; plot(B2023(11).neu_iu(:,1), B2023(11).neu_iu(:,4),'.');
%figure; plot(B2023(3).neu_iu(:,1), B2023(3).neu_iu(:,3),'.'); title('qiet - see 178')

%% ATTEMPT 2: remove datapoints with UP (Z-position) ERROR > 0.0425 [m]
% check figure: temporal structure of ENU errors
figure; clf; plot(A2023(10).neu(:,2),A2023(10).neu(:,8)); % Up 1-sigma
hold on; xlim([190 230])
plot(A2023(10).neu(:,2),A2023(10).neu(:,6)); % East 1-sigma
plot(A2023(10).neu(:,2),A2023(10).neu(:,4)); % North 1-sigma
plot([100 300],[0.0425 0.0425],'-')
ylim([0 0.1]); xlim([190 200])
title('SQ21: 1\sigma errors on positions in 2022')
xlabel('Day of Year, 2022  [ UTC ]'); ylabel('TRACK 1\sigma errors  [ m ]')
legend('Up 1\sigma','East 1\sigma','North 1\sigma','Threshold level for Up 1\sigma')
grid on
%print(gcf,'-dpng','-r300',sprintf('catalogue_data_quality/error_temporal_structure_SQ21_2022_190_200_260116.png')); 

threshold_UP = 0.0425; % 0.0425 [ m ] 
for i=1:1:num_sta
    if i~=6
    % Flag rows with Up errors (column 7) > threshold_UP 
    FLAG = B2023(i).neu_iu(:,7) > threshold_UP;
    % if it's not a day of a lake drainge = 1:
    NO_DRAIN = round(B2023(i).neu_iu(:,1)) ~= days_EXEMPT;
    NO_DRAIN_ANY = NO_DRAIN(:,1).*NO_DRAIN(:,2).*NO_DRAIN(:,3).* ...
        NO_DRAIN(:,4).*NO_DRAIN(:,5).*NO_DRAIN(:,6); % lake drainages are 0's
    % FLAG FOR BOTH
    FLAG_BOTH = logical(FLAG.*NO_DRAIN_ANY); % will leave 1's for both flags
    % remove rows with Up errors > some threshold
    B2023(i).neu_iu(FLAG_BOTH,:) = []; % 
    % save number of removed data points to calculate data rejected
    R2022(i).data_rejected_threshold_UP(1,1) = sum(FLAG_BOTH);
    
    elseif i==6 % SQ13 fell over by 2022/195
    % Flag rows with Up errors (column 7) > threshold_UP 
    FLAG = B2023(i).neu_iu(1:137100,7) > threshold_UP;
    % if it's not a day of a lake drainge = 1:
    NO_DRAIN = round(B2023(i).neu_iu(1:137100,1)) ~= days_EXEMPT;
    NO_DRAIN_ANY = NO_DRAIN(:,1).*NO_DRAIN(:,2).*NO_DRAIN(:,3).* ...
        NO_DRAIN(:,4).*NO_DRAIN(:,5).*NO_DRAIN(:,6); % lake drainages are 0's
    % FLAG FOR BOTH
    FLAG_BOTH = logical(FLAG.*NO_DRAIN_ANY); % will leave 1's for both flags
    % remove rows with Up errors > some threshold
    B2023(i).neu_iu(FLAG_BOTH,:) = []; % 
    % save number of removed data points to calculate data rejected
    R2022(i).data_rejected_threshold_UP(1,1) = sum(FLAG_BOTH);
    else 
    end
end


%% proportion data points not used due to bias flags, outliers, threshold_UP
for i=1:num_sta
    % starting number of data points for DOY 170–230
    data_reject_totals(i,1) = sum(R2022(i).data_rejected_BF(:,2)); % starting # of data points
    data_reject_totals(i,2) = sum(R2022(i).data_rejected_BF(:,3)); % rejected due to BF
    data_reject_totals(i,3) = R2022(i).data_rejected_outliers(1,2); % rejected due to outliers
    data_reject_totals(i,4) = R2022(i).data_rejected_threshold_UP(1,1); % rejected due to threshold_up
    % rejected proportion
    data_reject_totals(i,5) = data_reject_totals(i,2)/data_reject_totals(i,1);
    data_reject_totals(i,6) = data_reject_totals(i,3)/data_reject_totals(i,1);
    data_reject_totals(i,7) = data_reject_totals(i,4)/data_reject_totals(i,1);
end
for i=6 % SQ13 fell over by day 2022/200
    % starting number of data points for DOY 145–200
    data_reject_totals(i,1) = sum(R2022(i).data_rejected_BF(1:50,2)); % starting # of data points
    data_reject_totals(i,2) = sum(R2022(i).data_rejected_BF(1:50,3)); % rejected due to BF
    data_reject_totals(i,3) = R2022(i).data_rejected_outliers(1,2); % rejected due to outliers
    data_reject_totals(i,4) = R2022(i).data_rejected_threshold_UP(1,1); % rejected due to threshold_up
    % rejected proportion
    data_reject_totals(i,5) = data_reject_totals(i,2)/data_reject_totals(i,1);
    data_reject_totals(i,6) = data_reject_totals(i,3)/data_reject_totals(i,1);
    data_reject_totals(i,7) = data_reject_totals(i,4)/137101;
end

% command-line print out for Methods section reporting:
data_reject_averages_BF = mean(data_reject_totals(:,5))
data_reject_averages_outliers = mean(data_reject_totals(:,6))
data_reject_averages_threshold_UP = mean(data_reject_totals(:,7))
data_reject_averages_all = data_reject_averages_BF+data_reject_averages_threshold_UP+data_reject_averages_outliers

% check figure: does rejected data amount vary across array? (answer: no)
figure; clf; 
plot(1:1:num_sta,data_reject_totals(:,5),'s-'); hold on
plot(1:1:num_sta,data_reject_totals(:,6),'o-');
plot(1:1:num_sta,data_reject_totals(:,7),'^-');
lgd = legend('TRACK bias flag > 2','outlier (>3\sigma from running mean)','TRACK Up 1\sigma estimate > 0.0425 m','Location','East');
ylabel('Proportion of Data Points Rejected in 2022/150–230'); xlabel('Station ID');
set(gca,'xticklabel',upper(station_names),'xtick',1:1:num_sta); 
grid on; xlim([0 22])
% print figure
%print(gcf,'-dpng','-r300',sprintf('catalogue_data_quality/data_rejected_2022_150_230_260116.png')); 


%% GNARLY EDGES STATION-BY-STATION, visualized in EAST data
gnarly_MHIH = [160.1 161.1; 165.9 166.5; 173.9 173.4; 177.4 178.15; 179.95 180.1; ...
            183.4 184.0; 184.4 184.6; 196.95 197.05;...
            204.25 204.3; 209.26 209.31; 225.95 226.05; 254.95 255.1; ...
            262.95 263.05; 330.4 330.8];
for i=1:length(gnarly_MHIH(:,2)) % MHIH
    % Flag rows with time in gnarly bit
    FLAG1 = B2023(1).neu_iu(:,1) > gnarly_MHIH(i,1);
    FLAG2 = B2023(1).neu_iu(:,1) < gnarly_MHIH(i,2);  
    FLAG = logical(FLAG1.*FLAG2);
    % remove rows within this time window
    B2023(1).neu_iu(FLAG,:) = [];  % 
end
%figure; plot(B2023(1).neu_iu(:,1), B2023(1).neu_iu(:,3),'.');

gnarly_MLOW = [150.0 150.2; 151.7 152.2; 167.8 168.0; 183.325 183.7; ...
    184.05 184.3; 185.7 185.8; 203 203.4; 209 209.3; 213.95 214.35;...
    221.2 221.55; 221.7 221.96; 240.2 240.45];
for i=1:length(gnarly_MLOW(:,2)) % MLOW
    % Flag rows with time in gnarly bit
    FLAG1 = B2023(2).neu_iu(:,1) > gnarly_MLOW(i,1);
    FLAG2 = B2023(2).neu_iu(:,1) < gnarly_MLOW(i,2);  
    FLAG = logical(FLAG1.*FLAG2);
    % remove rows within this time window
    B2023(2).neu_iu(FLAG,:) = [];  % 
end
%figure; plot(B2023(2).neu_iu(:,1), B2023(2).neu_iu(:,3),'.');

gnarly_QIET = [189.25 190; 204.1 204.5; 223.25 223.85];
for i=1:length(gnarly_QIET(:,2)) % QIET 
    % Flag rows with time in gnarly bit
    FLAG1 = B2023(3).neu_iu(:,1) > gnarly_QIET (i,1);
    FLAG2 = B2023(3).neu_iu(:,1) < gnarly_QIET (i,2);  
    FLAG = logical(FLAG1.*FLAG2);
    % remove rows within this time window
    B2023(3).neu_iu(FLAG,:) = [];  % 
end
%figure; plot(B2023(3).neu_iu(:,1), B2023(3).neu_iu(:,3),'.');

gnarly_SQ11 = [183.1 183.55; 196.825 196.87; 197.0 197.25; 200.92 200.95;...
                204.2 204.3; 308.0 308.8];
for i=1:length(gnarly_SQ11(:,1)) % SQ11
    % Flag rows with time in gnarly bit
    FLAG1 = B2023(4).neu_iu(:,1) > gnarly_SQ11(i,1);
    FLAG2 = B2023(4).neu_iu(:,1) < gnarly_SQ11(i,2);  
    FLAG = logical(FLAG1.*FLAG2);
    % remove rows within this time window
    B2023(4).neu_iu(FLAG,:) = [];  % 
end
%figure; plot(B2023(4).neu_iu(:,1), B2023(4).neu_iu(:,3),'.');

gnarly_SQ12 = [223.65 223.75; 239.9 240.4];
for i=1:length(gnarly_SQ12(:,2)) % SQ12
    % Flag rows with time in gnarly bit
    FLAG1 = B2023(5).neu_iu(:,1) > gnarly_SQ12(i,1);
    FLAG2 = B2023(5).neu_iu(:,1) < gnarly_SQ12(i,2);  
    FLAG = logical(FLAG1.*FLAG2);
    % remove rows within this time window
    B2023(5).neu_iu(FLAG,:) = [];  % 
end
%figure; plot(B2023(5).neu_iu(:,1), B2023(5).neu_iu(:,3),'.');

gnarly_SQ13 = [177.6 178.825; 200.88 200.96; 268.0 268.3; 269.1 269.35];
for i=1:length(gnarly_SQ13(:,2)) % SQ13
    % Flag rows with time in gnarly bit
    FLAG1 = B2023(6).neu_iu(:,1) > gnarly_SQ13(i,1);
    FLAG2 = B2023(6).neu_iu(:,1) < gnarly_SQ13(i,2);  
    FLAG = logical(FLAG1.*FLAG2);
    % remove rows within this time window
    B2023(6).neu_iu(FLAG,:) = [];  % 
end
%figure; plot(B2023(6).neu_iu(:,1), B2023(6).neu_iu(:,3),'.');

gnarly_SQ14 = [150 150.2; 153.05 153.1; 183.4 183.6; 204.1 204.3; ...
            209 209.3; 223.7 224.1; 240 240.4];
for i=1:length(gnarly_SQ14(:,2)) % SQ14
    % Flag rows with time in gnarly bit
    FLAG1 = B2023(7).neu_iu(:,1) > gnarly_SQ14(i,1);
    FLAG2 = B2023(7).neu_iu(:,1) < gnarly_SQ14(i,2);  
    FLAG = logical(FLAG1.*FLAG2);
    % remove rows within this time window
    B2023(7).neu_iu(FLAG,:) = [];  % 
end
%figure; plot(B2023(7).neu_iu(:,1), B2023(7).neu_iu(:,3),'.');

gnarly_SQ15 = [140.4 140.8];
for i=1:length(gnarly_SQ15(:,2)) % SQ15
    % Flag rows with time in gnarly bit
    FLAG1 = B2023(8).neu_iu(:,1) > gnarly_SQ15(i,1);
    FLAG2 = B2023(8).neu_iu(:,1) < gnarly_SQ15(i,2);  
    FLAG = logical(FLAG1.*FLAG2);
    % remove rows within this time window
    B2023(8).neu_iu(FLAG,:) = [];  % 
end
%figure; plot(B2023(8).neu_iu(:,1), B2023(8).neu_iu(:,3),'.');

gnarly_SQ16 = [204.1 204.4; 209 209.34; 223.6 224.1; 253.2 253.8];
for i=1:length(gnarly_SQ16(:,2)) % SQ16
    % Flag rows with time in gnarly bit
    FLAG1 = B2023(9).neu_iu(:,1) > gnarly_SQ16(i,1);
    FLAG2 = B2023(9).neu_iu(:,1) < gnarly_SQ16(i,2);  
    FLAG = logical(FLAG1.*FLAG2);
    % remove rows within this time window
    B2023(9).neu_iu(FLAG,:) = [];  % 
end
%figure; plot(B2023(9).neu_iu(:,1), B2023(9).neu_iu(:,3),'.');

gnarly_SQ21 = [167.8 168.05; 169 169.6; 171.6 172.3; 173 173.1; 192.93 192.96;...
        223.05 223.8];
for i=1:length(gnarly_SQ21(:,2)) % SQ21
    % Flag rows with time in gnarly bit
    FLAG1 = B2023(10).neu_iu(:,1) > gnarly_SQ21(i,1);
    FLAG2 = B2023(10).neu_iu(:,1) < gnarly_SQ21(i,2);  
    FLAG = logical(FLAG1.*FLAG2);
    % remove rows within this time window
    B2023(10).neu_iu(FLAG,:) = [];  % 
end
%figure; plot(B2023(10).neu_iu(:,1), B2023(10).neu_iu(:,3),'.');

gnarly_SQ22 = [196.7 197.05; 219.65 220; 220.7 220.8; 223.05 224.15];
for i=1:length(gnarly_SQ22(:,2)) % SQ22
    for j=11
    % Flag rows with time in gnarly bit
    FLAG1 = B2023(j).neu_iu(:,1) > gnarly_SQ22(i,1);
    FLAG2 = B2023(j).neu_iu(:,1) < gnarly_SQ22(i,2);  
    FLAG = logical(FLAG1.*FLAG2);
    % remove rows within this time window
    B2023(j).neu_iu(FLAG,:) = [];  % 
    end
end
%figure; plot(B2023(11).neu_iu(:,1), B2023(11).neu_iu(:,3),'.');

gnarly_SQ23 = [201.9 202.15; 219.7 219.8; 223.0 224.2];
for i=1:length(gnarly_SQ23(:,2)) % SQ23
    for j=12
    % Flag rows with time in gnarly bit
    FLAG1 = B2023(j).neu_iu(:,1) > gnarly_SQ23(i,1);
    FLAG2 = B2023(j).neu_iu(:,1) < gnarly_SQ23(i,2);  
    FLAG = logical(FLAG1.*FLAG2);
    % remove rows within this time window
    B2023(j).neu_iu(FLAG,:) = [];  % 
    end
end
%figure; plot(B2023(12).neu_iu(:,1), B2023(12).neu_iu(:,3),'.');

gnarly_SQ24 = [196.95 197.1; 223.0 224.2];
for i=1:length(gnarly_SQ24(:,2)) % SQ24
    for j=13
    % Flag rows with time in gnarly bit
    FLAG1 = B2023(j).neu_iu(:,1) > gnarly_SQ24(i,1);
    FLAG2 = B2023(j).neu_iu(:,1) < gnarly_SQ24(i,2);  
    FLAG = logical(FLAG1.*FLAG2);
    % remove rows within this time window
    B2023(j).neu_iu(FLAG,:) = [];  % 
    end
end
%figure; plot(B2023(13).neu_iu(:,1), B2023(13).neu_iu(:,3),'.');

gnarly_SQ25 = [189.6 190.1; 196.8 197.08; 210.65 211.1; 223.42 223.55; 242.1 242.55];
for i=1:length(gnarly_SQ25(:,2)) % SQ25
    % Flag rows with time in gnarly bit
    FLAG1 = B2023(14).neu_iu(:,1) > gnarly_SQ25(i,1);
    FLAG2 = B2023(14).neu_iu(:,1) < gnarly_SQ25(i,2);  
    FLAG = logical(FLAG1.*FLAG2);
    % remove rows within this time window
    B2023(14).neu_iu(FLAG,:) = [];  % 
end
%figure; plot(B2023(14).neu_iu(:,1), B2023(14).neu_iu(:,3),'.');

gnarly_SQ26 = [196.75 197.1; 219.7 219.85; 223.05 223.6];
for i=1:length(gnarly_SQ26(:,2)) % SQ26
    for j=15
    % Flag rows with time in gnarly bit
    FLAG1 = B2023(j).neu_iu(:,1) > gnarly_SQ26(i,1);
    FLAG2 = B2023(j).neu_iu(:,1) < gnarly_SQ26(i,2);  
    FLAG = logical(FLAG1.*FLAG2);
    % remove rows within this time window
    B2023(j).neu_iu(FLAG,:) = [];  % 
    end
end
%figure; plot(B2023(15).neu_iu(:,1), B2023(15).neu_iu(:,3),'.');

gnarly_SQ27 = [182.9 183.2; 223.05 224.2];
for i=1:length(gnarly_SQ27(:,2)) % SQ27
    for j=16
    % Flag rows with time in gnarly bit
    FLAG1 = B2023(j).neu_iu(:,1) > gnarly_SQ27(i,1);
    FLAG2 = B2023(j).neu_iu(:,1) < gnarly_SQ27(i,2);  
    FLAG = logical(FLAG1.*FLAG2);
    % remove rows within this time window
    B2023(j).neu_iu(FLAG,:) = [];  % 
    end
end
%figure; plot(B2023(16).neu_iu(:,1), B2023(16).neu_iu(:,3),'.');

gnarly_SQ31 = [176.6 177; 183.6 183.75; 185.8 185.9; 220.8 221.7];
for i=1:length(gnarly_SQ31(:,2)) % SQ31
    for j=17
    % Flag rows with time in gnarly bit
    FLAG1 = B2023(j).neu_iu(:,1) > gnarly_SQ31(i,1);
    FLAG2 = B2023(j).neu_iu(:,1) < gnarly_SQ31(i,2);  
    FLAG = logical(FLAG1.*FLAG2);
    % remove rows within this time window
    B2023(j).neu_iu(FLAG,:) = [];  % 
    end
end
%figure; plot(B2023(17).neu_iu(:,1), B2023(17).neu_iu(:,3),'.');

gnarly_SQ34 = [210.05 210.1; 219.66 219.78; 220.3 220.5];
for i=1:length(gnarly_SQ34(:,2)) % SQ34
    for j=18
    % Flag rows with time in gnarly bit
    FLAG1 = B2023(j).neu_iu(:,1) > gnarly_SQ34(i,1);
    FLAG2 = B2023(j).neu_iu(:,1) < gnarly_SQ34(i,2);  
    FLAG = logical(FLAG1.*FLAG2);
    % remove rows within this time window
    B2023(j).neu_iu(FLAG,:) = [];  % 
    end
end
%figure; plot(B2023(18).neu_iu(:,1), B2023(18).neu_iu(:,3),'.');

gnarly_SQ35 = [203.8 204.45; 207.9 208.15; 209.0 209.35];
for i=1:length(gnarly_SQ35(:,2)) % SQ35
    for j=19
    % Flag rows with time in gnarly bit
    FLAG1 = B2023(j).neu_iu(:,1) > gnarly_SQ35(i,1);
    FLAG2 = B2023(j).neu_iu(:,1) < gnarly_SQ35(i,2);  
    FLAG = logical(FLAG1.*FLAG2);
    % remove rows within this time window
    B2023(j).neu_iu(FLAG,:) = [];  % 
    end
end
%figure; plot(B2023(19).neu_iu(:,1), B2023(19).neu_iu(:,3),'.');

gnarly_SQ36 = [183.6 184.1; 220.2 220.5; 221.65 221.9];
for i=1:length(gnarly_SQ36(:,2)) % SQ35
    for j=20
    % Flag rows with time in gnarly bit
    FLAG1 = B2023(j).neu_iu(:,1) > gnarly_SQ36(i,1);
    FLAG2 = B2023(j).neu_iu(:,1) < gnarly_SQ36(i,2);  
    FLAG = logical(FLAG1.*FLAG2);
    % remove rows within this time window
    B2023(j).neu_iu(FLAG,:) = [];  % 
    end
end
%figure; plot(B2023(20).neu_iu(:,1), B2023(20).neu_iu(:,3),'.');

gnarly_SQ37 = [177.2 177.6; 180.7 181.6; 183.3 183.38; 223.42 223.8];
for i=1:length(gnarly_SQ37(:,2)) % SQ35
    for j=21
    % Flag rows with time in gnarly bit
    FLAG1 = B2023(j).neu_iu(:,1) > gnarly_SQ37(i,1);
    FLAG2 = B2023(j).neu_iu(:,1) < gnarly_SQ37(i,2);  
    FLAG = logical(FLAG1.*FLAG2);
    % remove rows within this time window
    B2023(j).neu_iu(FLAG,:) = [];  % 
    end
end
%figure; plot(B2023(21).neu_iu(:,1), B2023(21).neu_iu(:,3),'.');

%% save B structure of cleaned positions
B2022=B2023; 
save B2022R_BF2_UP4_sZEROm_OUT6hr_clean_260113.mat B2022
num_sta = length([B2022.site])

%% calculate velocities and displacements
% ************** SWITCH FROM N-E-U TO E-N-U ORDERING *********************
% time_evaluate = 0.05; % evaluate every X time [ days] (72 minutes)
time_evaluate = 30/(24*60); % evaluate every X time [ days] (30 minutes)
time_vec_daily = time_start:time_evaluate:time_end; % use the values at midday (less noisy than day boundaries) [ DOY ]

% sliding least-squares window width
LS_sliding_width = 18/24; % evaluate over X time window [ days ] --- 18 or 36 hours
midpoint = LS_sliding_width./2; % [ days ]

% Data availability thresholds
% "have at least X hours of data points to calculate a linear fit to that segment of data"
THRESHOLD = (4*60)*6; % [ points ] (240 data points in an hour)*(hours that you want) --- 6 or 12 hours

% NOTES on CONSTRAINTS
% 'loose' = 6-hr remove outliers; or 0-hr smoothing;  %% STANDARD
%           36-hr sliding LS window; >=12hr of data points.
%
% 'tight' = 6-hr remove outliers; or 0-hr smoothing; %% DECENT!
%           18-hr sliding LS window; >=6hr of data points.

clear daily 

for j=1:1:length(time_vec_daily)-1 
    for i=1:1:num_sta % 100s, 200s, 300s
    
    % useful time information to keep
    daily(i,1).time_evaluate = time_evaluate; % evaluate every X time [ days]
    daily(i,1).LS_sliding_width = LS_sliding_width; % evaluate over X time window [ decimal days ]
    daily(i,1).t22(j,1) = time_vec_daily(j)+midpoint; % time
        
    % time indices (now pre-allocating)
    idx_time(i,j,1) = find((time_vec_daily(j)-B2023(i).neu_iu(:,1))<=0,1); % start
    idx_time(i,j,2) = find((time_vec_daily(j)+LS_sliding_width-B2023(i).neu_iu(:,1))<=0,1); % stop

    % displacement [ linear regression ]
    available_data = length(B2023(i).neu_iu(idx_time(i,j,1):idx_time(i,j,2),1));
        
        if available_data >= THRESHOLD % have at least THRESHOLD hours of data points to calculate a linear fit for the 'day'
        % slope
        [s_east, Err] = polyfit(B2023(i).neu_iu(idx_time(i,j,1):idx_time(i,j,2),1),...
            B2023(i).neu_iu(idx_time(i,j,1):idx_time(i,j,2),3),1); % slope in east [m/day]
        east_err = sqrt(diag((Err.R \ inv(Err.R')) / Err.normr^2 / Err.df));
        s_east_err = east_err(1); % slope error in east [m/day]
        [s_north, Err] = polyfit(B2023(i).neu_iu(idx_time(i,j,1):idx_time(i,j,2),1),...
            B2023(i).neu_iu(idx_time(i,j,1):idx_time(i,j,2),2),1); % slope in north [m/day]
        north_err = sqrt(diag((Err.R \ inv(Err.R')) / Err.normr^2 / Err.df));
        s_north_err = north_err(1); % slope error in north [m/day]
        [s_up, Err] = polyfit(B2023(i).neu_iu(idx_time(i,j,1):idx_time(i,j,2),1),...
            B2023(i).neu_iu(idx_time(i,j,1):idx_time(i,j,2),4),1); % slope in up [m/day]
        up_err = sqrt(diag((Err.R \ inv(Err.R')) / Err.normr^2 / Err.df));
        s_up_err = up_err(1); % slope error in north [m/day]

        % displacement from slope
        daily(i,1).displ_lin_fit(j,1) = s_east(1); % E [ m/d ]
        daily(i,1).displ_lin_fit(j,2) = s_north(1); % N [ m/d ]
        daily(i,1).displ_lin_fit(j,3) = s_up(1); % U [ m/d ] (z_s)
        daily(i,1).displ_lin_fit(j,4) = sqrt((daily(i,1).displ_lin_fit(j,1).^2)+(daily(i,1).displ_lin_fit(j,2).^2)); % horizontal (u_s) [ m/d ]
        
        % velocity from slope
        daily(i,1).vel_lin_fit(j,1) = s_east(1)*365.25; % [ m/yr ] E
        daily(i,1).vel_lin_fit(j,2) = s_north(1)*365.25; % [ m/yr ] N
        daily(i,1).vel_lin_fit(j,3) = s_up(1)*365.25; % [ m/yr ] U vertical (w_s)
        daily(i,1).vel_lin_fit(j,4) = sqrt((daily(i,1).displ_lin_fit(j,1).^2)+(daily(i,1).displ_lin_fit(j,2).^2)).*365.25; %  horizontal (u_s) [ m/yr ]
        
        % error on velocity estimate 
        daily(i,1).vel_lin_fit(j,5) = s_east_err*365.25; % [ m/yr ] E error
        daily(i,1).vel_lin_fit(j,6) = s_north_err*365.25; % [ m/yr ] N error
        daily(i,1).vel_lin_fit(j,7) = s_up_err*365.25; % [ m/yr ] U vertical error (w_s)

        % position [ from LS_sliding_width window average centred on the point ] in Euclidean space
        % position at time_vec_daily(1) = the starting point in time.
            if i<=17 
            daily(i,1).pos_24hrs(j,1) = (environs.X_km(i,1)*1e3) + ...
                nanmean(B2023(i).neu_iu(idx_time(i,j,1):idx_time(i,j,2),3)) - B2023(i).neu_iu(idx_time(i,1,1),3); % [ m ] E position
            daily(i,1).pos_24hrs(j,2) = (environs.Y_km(i,1)*1e3) + ...
                nanmean(B2023(i).neu_iu(idx_time(i,j,1):idx_time(i,j,2),2)) - B2023(i).neu_iu(idx_time(i,1,1),2); % [ m ] N position
            daily(i,1).pos_24hrs(j,3) = nanmean(B2023(i).neu_iu(idx_time(i,j,1):idx_time(i,j,2),4)) - B2023(i).neu_iu(idx_time(i,1,1),4); % [ m ] U position
            else % skip over SQ33 (#18 in environs) ONLY (SQ32 not in environs.mat in 2022)
            daily(i,1).pos_24hrs(j,1) = (environs.X_km(i+1,1)*1e3) + ...
                nanmean(B2023(i).neu_iu(idx_time(i,j,1):idx_time(i,j,2),3)) - B2023(i).neu_iu(idx_time(i,1,1),3); % [ m ] E position
            daily(i,1).pos_24hrs(j,2) = (environs.Y_km(i+1,1)*1e3) + ...
                nanmean(B2023(i).neu_iu(idx_time(i,j,1):idx_time(i,j,2),2)) - B2023(i).neu_iu(idx_time(i,1,1),2); % [ m ] N position
            daily(i,1).pos_24hrs(j,3) = nanmean(B2023(i).neu_iu(idx_time(i,j,1):idx_time(i,j,2),4)) - B2023(i).neu_iu(idx_time(i,1,1),4); % [ m ] U position    
            end    
                
        else % throw NaN when <THRESHOLD hours of data points in the time interval
        
        % displacement from slope
        daily(i,1).displ_lin_fit(j,1:4) = NaN; % [ m/interval ]
        % velocity from slope
        daily(i,1).vel_lin_fit(j,1:7) = NaN; % [ m/yr ] 
        % position [ from 24-hr average centred on the point ] in Euclidean space
        % position at time_vec_daily(1) = the starting point in time.
        daily(i,1).pos_24hrs(j,1:3) = NaN; 
        
        end
    end
end

% calculate strain rates for 100s, 200s, 300s and tie points
for i=1:1:num_sta % 100s, 200s, 300s
    for j=1:1:num_sta % 100s, 200s, 300s
        for t=1:1:length(time_vec_daily)-1
        % difference velocity of stations
        strain_rates.xvel_diff(i,j,t) = daily(i,1).vel_lin_fit(t,1)-daily(j,1).vel_lin_fit(t,1); % [ m/yr ] E  %%%% 2024-03-21 ERROR WAS HERE %%%%%
        strain_rates.yvel_diff(i,j,t) = daily(i,1).vel_lin_fit(t,2)-daily(j,1).vel_lin_fit(t,2); % [ m/yr ] N  %%%% 2024-03-21 ERROR WAS HERE %%%%%
        strain_rates.zvel_diff(i,j,t) = daily(i,1).vel_lin_fit(t,3)-daily(j,1).vel_lin_fit(t,3); % [ m/yr ] U  %%%% 2024-03-21 ERROR WAS HERE %%%%%
        % difference position of stations
        strain_rates.xpos_diff(i,j,t) = daily(i,1).pos_24hrs(t,1)-daily(j,1).pos_24hrs(t,1); % [ m ] E position
        strain_rates.ypos_diff(i,j,t) = daily(i,1).pos_24hrs(t,2)-daily(j,1).pos_24hrs(t,2); % [ m ] N position
        strain_rates.zpos_diff(i,j,t) = daily(i,1).pos_24hrs(t,3)-daily(j,1).pos_24hrs(t,3); % [ m ] U position;
        end
    end
end

% upper triangle: select only upper triangle of matrix U = triu(X)
lengthpos = length(strain_rates.xpos_diff);
for i=1:1:lengthpos
    % upper triangle of velocites
    du22(:,:,i)=triu(strain_rates.xvel_diff(:,:,i));
    dv22(:,:,i)=triu(strain_rates.yvel_diff(:,:,i));
    dw22(:,:,i)=triu(strain_rates.zvel_diff(:,:,i));
    % upper triangle of positions
    dx22(:,:,i)=triu(strain_rates.xpos_diff(:,:,i));
    dy22(:,:,i)=triu(strain_rates.ypos_diff(:,:,i));
    dz22(:,:,i)=triu(strain_rates.zpos_diff(:,:,i));
end

% 100s: Get positions and velocities into the along- and across-flow
% directions, when directions defined by 100's stations background
flowline_100s = nanmean(environs.u_s_bg_flow_dir_deg(4:9)) - 180; % 277.4585 --> 97.4585
dx22_f = dx22.*sind(flowline_100s) + dy22.*cosd(flowline_100s);
dy22_f = dy22.*sind(flowline_100s) + dx22.*cosd(flowline_100s);
du22_f = du22.*sind(flowline_100s) + dv22.*cosd(flowline_100s);
dv22_f = dv22.*sind(flowline_100s) + du22.*cosd(flowline_100s);

% STRAIN RATES
% calculate strain rates (velocity over length) [ year^{-1} ]
dudx_yr = du22./dx22; % east-west (~longitudinal) [ year^{-1} ]
dvdy_yr = dv22./dy22; % north-south (~transverse) [ year^{-1} ]
dwdz_yr = dw22./dz22; % vertical [ year^{-1} ]
shear_approx_yr = 0.5.*(((dv22)./dx22)+((du22)./dy22)); % shear [ year^{-1} ]
length_triu = length(dudx_yr);

% calculate strain rates (velocity over length) [year^{-1}] in direction of
% ice flow
lon_yr_100s = (du22_f)./dx22_f; % longitudinal [ year^{-1} ]
trans_yr_100s = (dv22_f)./dy22_f; % transverse [ year^{-1} ]
shear_yr_100s = 0.5.*(((dv22_f)./dx22_f)+((du22_f)./dy22_f)); % shear [ year^{-1} ]

% % ERRORS IN STRAIN RATES FOR 100s -- BROUGHT OVER FROM JGR:ES
% delta = 0.03; % assume all Track errors in horizontal = 3 cm (this is an overestimate in error of +1 cm)
% delta_sqr = delta*delta;
% this is all in units of year^{-1} 
for i=1:num_sta
    for j=1:num_sta
        for t=1:length(time_vec_daily)-1
           delta_lon_yr_100s(i,j,t) = abs(lon_yr_100s(i,j,t)).*sqrt(((1./(du22_f(i,j,t).^2)).*... % first half of error in velocity difference
               ((daily(i,1).vel_lin_fit(t,5)^2)+...  % errors in first station's E velocity
               (daily(j,1).vel_lin_fit(t,5)^2)))+... % errors in second station's E velocity
               ((1./(dx22_f(i,j,t).^2)).*((B2023(i).neu_iu(t,6)^2)+(B2023(j).neu_iu(t,6)^2)))); % final line is error in position difference 
           delta_trans_yr_100s(i,j,t) = abs(trans_yr_100s(i,j,t)).*sqrt(((1./(dv22_f(i,j,t).^2)).*... % first half of error in velocity difference
               ((daily(i,1).vel_lin_fit(t,6)^2)+...  % errors in first station's N velocity
               (daily(j,1).vel_lin_fit(t,6)^2)))+... % errors in second station's N velocity
               ((1./(dy22_f(i,j,t).^2)).*((B2023(i).neu_iu(t,5)^2)+(B2023(j).neu_iu(t,5)^2)))); % final line is error in position difference 
        end
    end
end

% 200's positions and velocities into the along- and across-flow 
% directions, when directions defined by 200's stations background
flowline_200s = nanmean(environs.u_s_bg_flow_dir_deg(10:16)) - 180; % 85.9536
dx22_f = dx22.*sind(flowline_200s) + dy22.*cosd(flowline_200s);
dy22_f = dy22.*sind(flowline_200s) + dx22.*cosd(flowline_200s);
du22_f = du22.*sind(flowline_200s) + dv22.*cosd(flowline_200s);
dv22_f = dv22.*sind(flowline_200s) + du22.*cosd(flowline_200s);

% STRAIN RATES
% calculate strain rates (velocity over length) [ year^{-1} ]
dudx_yr = du22./dx22; % east-west (~longitudinal) [ year^{-1} ]
dvdy_yr = dv22./dy22; % north-south (~transverse) [ year^{-1} ]
dwdz_yr = dw22./dz22; % vertical [ year^{-1} ]
shear_approx_yr = 0.5.*(((dv22)./dx22)+((du22)./dy22)); % shear [ year^{-1} ]
length_triu = length(dudx_yr);

% calculate strain rates (velocity over length) [year^{-1}] in direction of
% ice flow
lon_yr_200s = (du22_f)./dx22_f; % longitudinal [ year^{-1} ]
trans_yr_200s = (dv22_f)./dy22_f; % transverse [ year^{-1} ]
shear_yr_200s = 0.5.*(((dv22_f)./dx22_f)+((du22_f)./dy22_f)); % shear [ year^{-1} ]

% % ERRORS IN STRAIN RATES FOR 200s -- BROUGHT OVER FROM JGR:ES
% this is all in units of year^{-1} 
for i=1:num_sta
    for j=1:num_sta
        for t=1:length(time_vec_daily)-1
           % delta_lon)            
           delta_lon_yr_200s(i,j,t) = abs(lon_yr_200s(i,j,t)).*sqrt(((1./(du22_f(i,j,t).^2)).*... % first half of error in velocity difference
               ((daily(i,1).vel_lin_fit(t,5)^2)+...  % errors in first station's E velocity
               (daily(j,1).vel_lin_fit(t,5)^2)))+... % errors in second station's E velocity
               ((1./(dx22_f(i,j,t).^2)).*((B2023(i).neu_iu(t,6)^2)+(B2023(j).neu_iu(t,6)^2)))); % final line is error in position difference 
          delta_trans_yr_200s(i,j,t) = abs(trans_yr_200s(i,j,t)).*sqrt(((1./(dv22_f(i,j,t).^2)).*... % first half of error in velocity difference
               ((daily(i,1).vel_lin_fit(t,6)^2)+...  % errors in first station's N velocity
               (daily(j,1).vel_lin_fit(t,6)^2)))+... % errors in second station's N velocity
               ((1./(dy22_f(i,j,t).^2)).*((B2023(i).neu_iu(t,5)^2)+(B2023(j).neu_iu(t,5)^2)))); % final line is error in position difference 
        end
    end
end

% 300's positions and velocities into the along- and across-flow 
% directions, when directions defined by 300's stations background
flowline_300s = nanmean(environs.u_s_bg_flow_dir_deg(17:end)) - 180; % 97.0406
dx22_f = dx22.*sind(flowline_300s) + dy22.*cosd(flowline_300s);
dy22_f = dy22.*sind(flowline_300s) + dx22.*cosd(flowline_300s);
du22_f = du22.*sind(flowline_300s) + dv22.*cosd(flowline_300s);
dv22_f = dv22.*sind(flowline_300s) + du22.*cosd(flowline_300s);

% STRAIN RATES
% calculate strain rates (velocity over length) [ year^{-1} ]
dudx_yr = du22./dx22; % east-west (~longitudinal) [ year^{-1} ]
dvdy_yr = dv22./dy22; % north-south (~transverse) [ year^{-1} ]
dwdz_yr = dw22./dz22; % vertical [ year^{-1} ]
shear_approx_yr = 0.5.*(((dv22)./dx22)+((du22)./dy22)); % shear [ year^{-1} ]
length_triu = length(dudx_yr);

% calculate strain rates (velocity over length) [year^{-1}] in direction of
% ice flow
lon_yr_300s = (du22_f)./dx22_f; % longitudinal [ year^{-1} ]
trans_yr_300s = (dv22_f)./dy22_f; % transverse [ year^{-1} ]
shear_yr_300s = 0.5.*(((dv22_f)./dx22_f)+((du22_f)./dy22_f)); % shear [ year^{-1} ]

% % ERRORS IN STRAIN RATES FOR 300s -- BROUGHT OVER FROM JGR:ES
% this is all in units of year^{-1} 
for i=1:num_sta
    for j=1:num_sta
        for t=1:length(time_vec_daily)-1
           % delta_lon)            
           % delta_lon_yr_300s(i,j,t) = abs(lon_yr_300s(i,j,t)).*sqrt(((1./(du22_f(i,j,t).^2)).*... % first half of error in velocity difference
           %     ((2*daily(i,1).vel_lin_fit(t,5))+(2*daily(j,1).vel_lin_fit(t,5))))+... % errors in both stations' E velocities
           %     ((1./(dx22_f(i,j,t).^2)).*(2*delta_sqr))); % final line is error in displacement
            delta_lon_yr_300s(i,j,t) = abs(lon_yr_300s(i,j,t)).*sqrt(((1./(du22_f(i,j,t).^2)).*... % first half of error in velocity difference
               ((daily(i,1).vel_lin_fit(t,5)^2)+...  % errors in first station's E velocity
               (daily(j,1).vel_lin_fit(t,5)^2)))+... % errors in second station's E velocity
               ((1./(dx22_f(i,j,t).^2)).*((B2023(i).neu_iu(t,6)^2)+(B2023(j).neu_iu(t,6)^2)))); % final line is error in position difference 
           delta_trans_yr_300s(i,j,t) = abs(trans_yr_300s(i,j,t)).*sqrt(((1./(dv22_f(i,j,t).^2)).*... % first half of error in velocity difference
               ((daily(i,1).vel_lin_fit(t,6)^2)+...  % errors in first station's N velocity
               (daily(j,1).vel_lin_fit(t,6)^2)))+... % errors in second station's N velocity
               ((1./(dy22_f(i,j,t).^2)).*((B2023(i).neu_iu(t,5)^2)+(B2023(j).neu_iu(t,5)^2)))); % final line is error in position difference 
        end
    end
end

%% strain rates save on their own
daily_strain_rates_2023.station_names = station_names'; % station names
daily_strain_rates_2023.time = daily(1).t22; % time
daily_strain_rates_2023.lon_yr_100s = lon_yr_100s;
daily_strain_rates_2023.lon_yr_200s = lon_yr_200s;
daily_strain_rates_2023.lon_yr_300s = lon_yr_300s;
daily_strain_rates_2023.trans_yr_100s = trans_yr_100s;
daily_strain_rates_2023.trans_yr_200s = trans_yr_200s;
daily_strain_rates_2023.trans_yr_300s = trans_yr_300s;
daily_strain_rates_2023.shear_yr_100s = shear_yr_100s;
daily_strain_rates_2023.shear_yr_200s = shear_yr_200s;
daily_strain_rates_2023.shear_yr_300s = shear_yr_300s;
% also save the errors
daily_strain_rates_2023.delta_lon_yr_100s = delta_lon_yr_100s;
daily_strain_rates_2023.delta_lon_yr_200s = delta_lon_yr_200s;
daily_strain_rates_2023.delta_lon_yr_300s = delta_lon_yr_300s;
daily_strain_rates_2023.delta_trans_yr_100s = delta_trans_yr_100s;
daily_strain_rates_2023.delta_trans_yr_200s = delta_trans_yr_200s;
daily_strain_rates_2023.delta_trans_yr_300s = delta_trans_yr_300s;
% slideing least squares velocities & their errors (E,N,U,vel,delta-E,delta-N,delta-U)
for i=1:num_sta
daily_strain_rates_2023.vel_lin_fit(i,1:length(daily(1).t22),1:7) = daily(i).vel_lin_fit; 
end

%save daily_strain_rates_2022R_30min_BF2_UP4_sZERO_clean_w18_t6_260119.mat daily_strain_rates_2023

%% daily c_dot_delta_t = [vertical displacement] less [advective] less [strain]
% calculate epsilon_dot_zz from station velocities, pick stations to pair
% need to make a figure of these choices of station pairs -- it's
% subjective based on which pairs are least likely to be swamped by
% lake-drainage events

% % 1 = MHIH; 2 = MLOW; 3=QIET; 
% % 4=SQ11;   5=SQ12;   6=SQ13;   7=SQ14;   8=SQ15;   9=SQ16;
% % 10=SQ21;  11=SQ22;  12=SQ23;  13=SQ24;  14=SQ25;  15=SQ26;  16=SQ27;
% % 17=SQ31;  18=SQ34;  19=SQ35;  20=SQ36;  21=SQ37;

% TIES
% MHIH
% lon = MHIH-SQ31 (1,17); trans = MHIH-SQ22 (1,11)
epsilon_dot(1,:).epsilon_dot_lon = squeeze(lon_yr_200s(1,17,:)); % [ year^{-1} ]
epsilon_dot(1,:).epsilon_dot_trans = squeeze(trans_yr_200s(1,11,:)); % [ year^{-1} ]
epsilon_dot(1,:).epsilon_dot_zz = -1.*(epsilon_dot(1,:).epsilon_dot_lon+epsilon_dot(1,:).epsilon_dot_trans); % [ year^{-1} ]
epsilon_dot(1,:).epsilon_dot_zz_err = (squeeze(delta_lon_yr_200s(1,17,:)).^2)+(squeeze(delta_trans_yr_200s(1,11,:)).^2); % [ year^{-1} ]
% MLOW
% lon = MLOW-QIET (2,3); trans = MLOW-QIET (2,3)
epsilon_dot(2,:).epsilon_dot_lon = squeeze(lon_yr_200s(2,14,:)); % [ year^{-1} ]
epsilon_dot(2,:).epsilon_dot_trans = squeeze(trans_yr_200s(2,3,:)); % [ year^{-1} ]
epsilon_dot(2,:).epsilon_dot_zz = -1.*(epsilon_dot(2,:).epsilon_dot_lon+epsilon_dot(2,:).epsilon_dot_trans); % [ year^{-1} ]
epsilon_dot(2,:).epsilon_dot_zz_err = (squeeze(delta_lon_yr_200s(2,14,:)).^2)+(squeeze(delta_trans_yr_200s(2,3,:)).^2); % [ year^{-1} ]
% QIET
% lon = QIET-SQ12 (3,5); trans = QIET-SQ16 (3,9)
epsilon_dot(3,:).epsilon_dot_lon = squeeze(lon_yr_200s(3,5,:)); % [ year^{-1} ]
epsilon_dot(3,:).epsilon_dot_trans = squeeze(trans_yr_200s(3,9,:)); % [ year^{-1} ]
epsilon_dot(3,:).epsilon_dot_zz = -1.*(epsilon_dot(3,:).epsilon_dot_lon+epsilon_dot(3,:).epsilon_dot_trans); % [ year^{-1} ]
epsilon_dot(3,:).epsilon_dot_zz_err = (squeeze(delta_lon_yr_200s(3,5,:)).^2)+(squeeze(delta_trans_yr_200s(3,9,:)).^2); % [ year^{-1} ]


% 100s
% SQ11
% lon = SQ11-SQ12 (4,5); trans = SQ11-SQ14 (4,7)
epsilon_dot(4,:).epsilon_dot_lon = squeeze(lon_yr_100s(4,5,:)); % [ year^{-1} ]
epsilon_dot(4,:).epsilon_dot_trans = squeeze(trans_yr_100s(4,7,:)); % [ year^{-1} ]
epsilon_dot(4,:).epsilon_dot_zz = -1.*(epsilon_dot(4,:).epsilon_dot_lon+epsilon_dot(4,:).epsilon_dot_trans); % [ year^{-1} ]
epsilon_dot(4,:).epsilon_dot_zz_err = (squeeze(delta_lon_yr_100s(4,5,:)).^2)+(squeeze(delta_trans_yr_100s(4,7,:)).^2); % [ year^{-1} ]
% SQ12
% lon = QIET-SQ12 (3,5); trans = SQ12-SQ14 (5,7)
epsilon_dot(5,:).epsilon_dot_lon = squeeze(lon_yr_100s(3,5,:)); % [ year^{-1} ]
epsilon_dot(5,:).epsilon_dot_trans = squeeze(trans_yr_100s(5,7,:)); % [ year^{-1} ]
epsilon_dot(5,:).epsilon_dot_zz = -1.*(epsilon_dot(5,:).epsilon_dot_lon+epsilon_dot(5,:).epsilon_dot_trans); % [ year^{-1} ]
epsilon_dot(5,:).epsilon_dot_zz_err = (squeeze(delta_lon_yr_100s(3,5,:)).^2)+(squeeze(delta_trans_yr_100s(5,7,:)).^2); % [ year^{-1} ]
% SQ13
% lon = QIET-SQ13 (3,6); trans = SQ13-SQ16 (6,9)
epsilon_dot(6,:).epsilon_dot_lon = squeeze(lon_yr_100s(3,6,:)); % [ year^{-1} ]
epsilon_dot(6,:).epsilon_dot_trans = squeeze(trans_yr_100s(6,9,:)); % [ year^{-1} ]
epsilon_dot(6,:).epsilon_dot_zz = -1.*(epsilon_dot(6,:).epsilon_dot_lon+epsilon_dot(6,:).epsilon_dot_trans); % [ year^{-1} ]
epsilon_dot(6,:).epsilon_dot_zz_err = (squeeze(delta_lon_yr_100s(3,6,:)).^2)+(squeeze(delta_trans_yr_100s(6,9,:)).^2); % [ year^{-1} ]
% SQ14
% lon = SQ12-SQ14 (5,7); trans = SQ14-SQ15 (7,8)
epsilon_dot(7,:).epsilon_dot_lon = squeeze(lon_yr_100s(5,7,:)); % [ year^{-1} ]
epsilon_dot(7,:).epsilon_dot_trans = squeeze(trans_yr_100s(7,8,:)); % [ year^{-1} ]
epsilon_dot(7,:).epsilon_dot_zz = -1.*(epsilon_dot(7,:).epsilon_dot_lon+epsilon_dot(7,:).epsilon_dot_trans); % [ year^{-1} ]
epsilon_dot(7,:).epsilon_dot_zz_err = (squeeze(delta_lon_yr_100s(5,7,:)).^2)+(squeeze(delta_trans_yr_100s(7,8,:)).^2); % [ year^{-1} ]
% SQ15
% lon = SQ15-SQ16 (8,9); trans = SQ14-SQ15 (7,8)
epsilon_dot(8,:).epsilon_dot_lon = squeeze(lon_yr_100s(8,9,:)); % [ year^{-1} ]
epsilon_dot(8,:).epsilon_dot_trans = squeeze(trans_yr_100s(7,8,:)); % [ year^{-1} ]
epsilon_dot(8,:).epsilon_dot_zz = -1.*(epsilon_dot(8,:).epsilon_dot_lon+epsilon_dot(8,:).epsilon_dot_trans); % [ year^{-1} ]
epsilon_dot(8,:).epsilon_dot_zz_err = (squeeze(delta_lon_yr_100s(8,9,:)).^2)+(squeeze(delta_trans_yr_100s(7,8,:)).^2); % [ year^{-1} ]
% SQ16
% lon = QIET-SQ16 (3,9); trans = SQ14-SQ16 (7,9)
epsilon_dot(9,:).epsilon_dot_lon = squeeze(lon_yr_100s(3,9,:)); % [ year^{-1} ]
epsilon_dot(9,:).epsilon_dot_trans = squeeze(trans_yr_100s(7,9,:)); % [ year^{-1} ]
epsilon_dot(9,:).epsilon_dot_zz = -1.*(epsilon_dot(9,:).epsilon_dot_lon+epsilon_dot(9,:).epsilon_dot_trans); % [ year^{-1} ]
epsilon_dot(9,:).epsilon_dot_zz_err = (squeeze(delta_lon_yr_100s(3,9,:)).^2)+(squeeze(delta_trans_yr_100s(7,9,:)).^2); % [ year^{-1} ]


% 200s
% SQ21
% lon = SQ21-SQ24 (10,13); trans = SQ21-SQ22 (10,11)
epsilon_dot(10,:).epsilon_dot_lon = squeeze(lon_yr_200s(10,13,:)); % [ year^{-1} ]
epsilon_dot(10,:).epsilon_dot_trans = squeeze(trans_yr_200s(10,11,:)); % [ year^{-1} ]
epsilon_dot(10,:).epsilon_dot_zz = -1.*(epsilon_dot(10,:).epsilon_dot_lon+epsilon_dot(10,:).epsilon_dot_trans); % [ year^{-1} ]
epsilon_dot(10,:).epsilon_dot_zz_err = (squeeze(delta_lon_yr_200s(10,13,:)).^2)+(squeeze(delta_trans_yr_200s(10,11,:)).^2); % [ year^{-1} ]
% SQ22
% lon = SQ22_SQ23 (11,12); trans = SQ22-SQ24 (11,13)
epsilon_dot(11,:).epsilon_dot_lon = squeeze(lon_yr_200s(11,12,:)); % [ year^{-1} ]
epsilon_dot(11,:).epsilon_dot_trans = squeeze(trans_yr_200s(11,13,:)); % [ year^{-1} ]
epsilon_dot(11,:).epsilon_dot_zz = -1.*(epsilon_dot(11,:).epsilon_dot_lon+epsilon_dot(11,:).epsilon_dot_trans); % [ year^{-1} ]
epsilon_dot(11,:).epsilon_dot_zz_err = (squeeze(delta_lon_yr_200s(11,12,:)).^2)+(squeeze(delta_trans_yr_200s(11,13,:)).^2); % [ year^{-1} ]
% SQ23
% lon = SQ23-SQ31 (12,17); trans = SQ23-SQ24 (12,13)
epsilon_dot(12,:).epsilon_dot_lon = squeeze(lon_yr_200s(12,17,:)); % [ year^{-1} ]
epsilon_dot(12,:).epsilon_dot_trans = squeeze(trans_yr_200s(12,13,:)); % [ year^{-1} ]
epsilon_dot(12,:).epsilon_dot_zz = -1.*(epsilon_dot(12,:).epsilon_dot_lon+epsilon_dot(12,:).epsilon_dot_trans); % [ year^{-1} ]
epsilon_dot(12,:).epsilon_dot_zz_err = (squeeze(delta_lon_yr_200s(12,17,:)).^2)+(squeeze(delta_trans_yr_200s(12,13,:)).^2); % [ year^{-1} ]
% SQ24
% lon = SQ24-SQ27 (13,16); trans = SQ24-SQ26 (13,15)
epsilon_dot(13,:).epsilon_dot_lon = squeeze(lon_yr_200s(13,16,:)); % [ year^{-1} ]
epsilon_dot(13,:).epsilon_dot_trans = squeeze(trans_yr_200s(13,15,:)); % [ year^{-1} ]
epsilon_dot(13,:).epsilon_dot_zz = -1.*(epsilon_dot(13,:).epsilon_dot_lon+epsilon_dot(13,:).epsilon_dot_trans); % [ year^{-1} ]
epsilon_dot(13,:).epsilon_dot_zz_err = (squeeze(delta_lon_yr_200s(13,16,:)).^2)+(squeeze(delta_trans_yr_200s(13,15,:)).^2); % [ year^{-1} ]
% SQ25
% lon = SQ25-SQ26 (14,15); trans = QIET-SQ25 (3,14)
epsilon_dot(14,:).epsilon_dot_lon = squeeze(lon_yr_200s(14,15,:)); % [ year^{-1} ]
epsilon_dot(14,:).epsilon_dot_trans = squeeze(trans_yr_200s(3,14,:)); % [ year^{-1} ]
epsilon_dot(14,:).epsilon_dot_zz = -1.*(epsilon_dot(14,:).epsilon_dot_lon+epsilon_dot(14,:).epsilon_dot_trans); % [ year^{-1} ]
epsilon_dot(14,:).epsilon_dot_zz_err = (squeeze(delta_lon_yr_200s(14,15,:)).^2)+(squeeze(delta_trans_yr_200s(3,14,:)).^2); % [ year^{-1} ]
% SQ26
% lon = SQ25-SQ26 (14,15); trans = QIET-SQ26 (3,15)
epsilon_dot(15,:).epsilon_dot_lon = squeeze(lon_yr_200s(14,15,:)); % [ year^{-1} ]
epsilon_dot(15,:).epsilon_dot_trans = squeeze(trans_yr_200s(3,15,:)); % [ year^{-1} ]
epsilon_dot(15,:).epsilon_dot_zz = -1.*(epsilon_dot(15,:).epsilon_dot_lon+epsilon_dot(15,:).epsilon_dot_trans); % [ year^{-1} ]
epsilon_dot(15,:).epsilon_dot_zz_err = (squeeze(delta_lon_yr_200s(14,15,:)).^2)+(squeeze(delta_trans_yr_200s(3,15,:)).^2); % [ year^{-1} ]
% SQ27
% lon = SQ27-SQ35 (16,20); trans = SQ26-SQ27 (15,16)
epsilon_dot(16,:).epsilon_dot_lon = squeeze(lon_yr_200s(16,20,:)); % [ year^{-1} ]
epsilon_dot(16,:).epsilon_dot_trans = squeeze(trans_yr_200s(15,16,:)); % [ year^{-1} ]
epsilon_dot(16,:).epsilon_dot_zz = -1.*(epsilon_dot(16,:).epsilon_dot_lon+epsilon_dot(16,:).epsilon_dot_trans); % [ year^{-1} ]
epsilon_dot(16,:).epsilon_dot_zz_err = (squeeze(delta_lon_yr_200s(16,20,:)).^2)+(squeeze(delta_trans_yr_200s(15,16,:)).^2); % [ year^{-1} ]


% 300s
% SQ31
% lon = SQ31-SQ34 (17,18); trans = SQ31-SQ34 (17,18)
epsilon_dot(17,:).epsilon_dot_lon = squeeze(lon_yr_300s(17,18,:)); % [ year^{-1} ]
epsilon_dot(17,:).epsilon_dot_trans = squeeze(trans_yr_300s(17,18,:)); % [ year^{-1} ]
epsilon_dot(17,:).epsilon_dot_zz = -1.*(epsilon_dot(17,:).epsilon_dot_lon+epsilon_dot(17,:).epsilon_dot_trans); % [ year^{-1} ]
epsilon_dot(17,:).epsilon_dot_zz_err = (squeeze(delta_lon_yr_300s(17,18,:)).^2)+(squeeze(delta_trans_yr_300s(17,18,:)).^2); % [ year^{-1} ]
% SQ34
% lon = SQ34-SQ37 (18,21); trans = SQ34-SQ35 (18,19)
epsilon_dot(18,:).epsilon_dot_lon = squeeze(lon_yr_300s(18,21,:)); % [ year^{-1} ]
epsilon_dot(18,:).epsilon_dot_trans = squeeze(trans_yr_300s(18,19,:)); % [ year^{-1} ]
epsilon_dot(18,:).epsilon_dot_zz = -1.*(epsilon_dot(18,:).epsilon_dot_lon+epsilon_dot(18,:).epsilon_dot_trans); % [ year^{-1} ]
epsilon_dot(18,:).epsilon_dot_zz_err = (squeeze(delta_lon_yr_300s(18,21,:)).^2)+(squeeze(delta_trans_yr_300s(18,19,:)).^2); % [ year^{-1} ]
% SQ35
% lon = SQ35-SQ36 (19,20); trans = SQ35-SQ37 (19,21)
epsilon_dot(19,:).epsilon_dot_lon = squeeze(lon_yr_300s(19,20,:)); % [ year^{-1} ]
epsilon_dot(19,:).epsilon_dot_trans = squeeze(trans_yr_300s(19,21,:)); % [ year^{-1} ]
epsilon_dot(19,:).epsilon_dot_zz = -1.*(epsilon_dot(19,:).epsilon_dot_lon+epsilon_dot(19,:).epsilon_dot_trans); % [ year^{-1} ]
epsilon_dot(19,:).epsilon_dot_zz_err = (squeeze(delta_lon_yr_300s(19,20,:)).^2)+(squeeze(delta_trans_yr_300s(19,21,:)).^2); % [ year^{-1} ]
% SQ36
% lon = SQ35-SQ36 (19,20); trans = SQ34-SQ36 (18,20)
epsilon_dot(20,:).epsilon_dot_lon = squeeze(lon_yr_300s(19,20,:)); % [ year^{-1} ]
epsilon_dot(20,:).epsilon_dot_trans = squeeze(trans_yr_300s(18,20,:)); % [ year^{-1} ]
epsilon_dot(20,:).epsilon_dot_zz = -1.*(epsilon_dot(20,:).epsilon_dot_lon+epsilon_dot(20,:).epsilon_dot_trans); % [ year^{-1} ]
epsilon_dot(20,:).epsilon_dot_zz_err = (squeeze(delta_lon_yr_300s(19,20,:)).^2)+(squeeze(delta_trans_yr_300s(18,20,:)).^2); % [ year^{-1} ]
% SQ37
% lon = SQ34-SQ37 (18,21); trans = SQ36-SQ37 (20,21)
epsilon_dot(21,:).epsilon_dot_lon = squeeze(lon_yr_300s(18,21,:)); % [ year^{-1} ]
epsilon_dot(21,:).epsilon_dot_trans = squeeze(trans_yr_300s(20,21,:)); % [ year^{-1} ]
epsilon_dot(21,:).epsilon_dot_zz = -1.*(epsilon_dot(21,:).epsilon_dot_lon+epsilon_dot(21,:).epsilon_dot_trans); % [ year^{-1} ]
epsilon_dot(21,:).epsilon_dot_zz_err = (squeeze(delta_lon_yr_300s(18,21,:)).^2)+(squeeze(delta_trans_yr_300s(20,21,:)).^2); % [ year^{-1} ]

%% calculate theta_b with epsilon_dot_zz from 2023/150-160 
% positive theta_b here = bed is sloping uphill in direction of ice flow
zz_bg_start_stop = [150, 160]; 
for i=1:1:num_sta
    % time indices
    idx_time(i,1) = find((zz_bg_start_stop(1)-daily(i).t22)<=0,1); % start
    idx_time(i,2) = find((zz_bg_start_stop(2)-daily(i).t22)<=0,1); % stop
    
    epsilon_dot_zz_bg_station_pairs(i,1) = nanmean(epsilon_dot(i,:).epsilon_dot_zz(idx_time(i,1):idx_time(i,2),1)); % [1/yr]
    theta_b_from_epsilon_zz(i,1) = atand((w_s_bg(i)-(epsilon_dot_zz_bg_station_pairs(i,1)*H(i)))./(u_b_bg(i))); % [deg]
end

% daily c_dot_delta_t = [vertical displacement] less [advective] less [strain]
% using epsilon_zz from station pairs
% discrete integral time interval
delta_t = LS_sliding_width; % [ days ]
delta_t_years = delta_t./365.25; % [ years ]
delta_t_years1 = 1/365.25;
% c_dot for given time interval
for i=1:1:num_sta % 100s, 200s, 300s
    % time array
    daily_epsilon_zz(i,1).delta_t = delta_t; % [ days ]  
    daily_epsilon_zz(i,1).t22 = daily(i).t22; % time vector
    daily_epsilon_zz(i,1).epsilon_dot_zz_bg_station_pairs = epsilon_dot_zz_bg_station_pairs; % epsilon_{zz,bg}
    daily_epsilon_zz(i,1).theta_b_from_epsilon_zz = theta_b_from_epsilon_zz; % theta_b
    daily_epsilon_zz(i,1).u_s_bg = u_s_bg; % u_b_bg = u_s_bg
    daily_epsilon_zz(i,1).u_s = daily(i,1).vel_lin_fit(:,4); % horizontal (u_s) [ m/yr ]
    daily_epsilon_zz(i,1).w_s = daily(i,1).vel_lin_fit(:,3); % vertical (w_s) [ m/yr ]

    % bed opening length [ m in one day ] [ instantaneous ]
    daily_epsilon_zz(i,1).c_dot_delta_t = (daily(i,1).displ_lin_fit(:,3).*365.25.*delta_t_years1) - ... % vertical dislacement z_s [m/day] [day/yr] [yr]
                       ((daily(i,1).vel_lin_fit(:,4)).*(tand(theta_b_from_epsilon_zz(i))*delta_t_years1)) - ... % advective [m/yr] [yr] --> [m]
                       (epsilon_dot(i,:).epsilon_dot_zz.*H(i).*delta_t_years1); % strain    [1/yr] [m] [yr] --> [m]
    
    % error in bed opening length [ m in one day ] [ instantaneous ]
    daily_epsilon_zz(i,1).c_dot_delta_t_err = sqrt(((daily(i,1).vel_lin_fit(:,7).*delta_t_years1).^2) + ... % vertical velocity error z_s [m/yr] [1/yr] --> [m]
                       ((daily(i,1).vel_lin_fit(:,5).*delta_t_years1.*deg2rad(theta_b_from_epsilon_zz(i))).^2) + ... % advective error 1 [m/yr] [yr] --> [m]
                       ((daily(i,1).vel_lin_fit(:,4).*delta_t_years1.*deg2rad(theta_b_err(i))).^2) + ... % advective error 2 [m/yr] [yr] --> [m]
                       ((epsilon_dot(i,:).epsilon_dot_zz.*delta_t_years1.*H_err(i)).^2) + ...  % strain error 1  [1/yr] [m] [yr] --> [m] % bed opening rate [ m/day ]
                       ((H(i).^2).*(epsilon_dot(i,:).epsilon_dot_zz_err.*delta_t_years1))); % strain error 2  [1/yr] [m] [yr] --> [m] % bed opening rate [ m/day ]
    % bed opening cumulative [ m ] (time zero is time_vec_daily(1))
    %daily_epsilon_zz(i,1).c_dot_delta_t_cumulative = cumsum(daily(i,1).c_dot_delta_t.*time_evaluate,"omitnan"); % [m/day][day] --> [m]
    % c_dot_delta_t components
    daily_epsilon_zz(i,1).c_dot_delta_t_vert_displ = daily(i,1).displ_lin_fit(:,3).*365.25.*delta_t_years1; % vertical dislacement z_s [m/yr]
    daily_epsilon_zz(i,1).c_dot_delta_t_advective = ((daily(i,1).vel_lin_fit(:,4)./u_s_bg(i)).*(u_b_bg(i).*tand(theta_b_from_epsilon_zz(i)).*delta_t_years1)); % advective [m/yr] [yr] --> [m]
    daily_epsilon_zz(i,1).c_dot_delta_t_epsilon_zz = epsilon_dot(i,:).epsilon_dot_zz.*H(i).*delta_t_years1; % strain [1/yr] [m] [yr] --> [m]
    % bed opening cumulative components [ m ] (time zero is time_vec_daily(1))
    daily_epsilon_zz(i,1).c_dot_delta_t_vert_displ_cumulative = cumsum(daily_epsilon_zz(i,1).c_dot_delta_t_vert_displ.*time_evaluate,"omitnan"); % [ m ]
    daily_epsilon_zz(i,1).c_dot_delta_t_advective_cumulative = cumsum(daily_epsilon_zz(i,1).c_dot_delta_t_advective.*time_evaluate,"omitnan"); % [ m ]
    daily_epsilon_zz(i,1).c_dot_delta_t_epsilon_zz_cumulative = cumsum(daily_epsilon_zz(i,1).c_dot_delta_t_epsilon_zz.*time_evaluate,"omitnan"); % [ m ]

    % stain rates
    daily_epsilon_zz(i,1).epsilon_dot_lon = epsilon_dot(i,:).epsilon_dot_lon; % [ year^{-1} ]
    daily_epsilon_zz(i,1).epsilon_dot_trans = epsilon_dot(i,:).epsilon_dot_trans; % [ year^{-1} ]
    daily_epsilon_zz(i,1).epsilon_dot_zz = epsilon_dot(i,:).epsilon_dot_zz; % [ year^{-1} ]
end

%save daily_epsilon_zz_2022R_BF2_30min_BF2_UP4_sZERO_clean_w18_t6_260119.mat daily_epsilon_zz

%% Figure checks to view the data
%%% $\dot{\epsilon}_{zz,bg}$
% 
% figure(1); clf; 
% set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',1.2.*[5 2 22 12]);
% axe1 = axes('Position',[0.08 0.1 0.9 0.85],'Box','on','NextPlot','add','XTickLabels',[]);
% axes(axe1)
% plot(epsilon_dot_zz_bg,'ks','MarkerSize',8); hold on;
% plot(epsilon_dot_zz_bg_station_pairs,'ro','MarkerSize',8); 
% for i=1:1:num_sta; plot([i i],[epsilon_dot_zz_bg(i) epsilon_dot_zz_bg_station_pairs(i)],'-k','LineWidth',1.1); end
% set(gca,'xtick',1:1:21,'xticklabel',station_names); axe1.YAxis.Exponent = 0;
% legend('from BedMachine \theta_{b}','from station-pair $\dot{\epsilon}_{zz,bg}$','Interpreter','latex','Location','SouthEast');
% %for i=1:1:22; text(i+0.2,epsilon_dot_zz_bg(i,1),sprintf('%5.4f',epsilon_dot_zz_bg(i,1))); end
% grid on; xlim([0 22]); %ylim([-0.01 0.01])
% title('thickness-integrated vertical strain rate during 2023/150-160  [ yr^{-1} ]')
% ylabel('$\dot{\epsilon}_{zz,bg}$  [ yr$^{-1}$ ]','Interpreter','latex','FontSize',14)
% 
% % print figure
% % cd catalogue_c_dot/30days_48hrs/
% % print(gcf,'-dpng','-r500',sprintf('epsilonzz_2023_150_155_compare_20240121.png')); 
% % cd ../../
% 
% %% $\theta_{b}$
% 
% figure(2); clf
% set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',1.2.*[5 2 22 12]);
% axe1 = axes('Position',[0.08 0.1 0.9 0.85],'Box','on','NextPlot','add','XTickLabels',[]);
% axes(axe1)
% plot(-1.*theta_b,'ks','MarkerSize',8); hold on;
% plot(theta_b_from_epsilon_zz,'ro','MarkerSize',8); 
% for i=1:1:num_sta; plot([i i],[-1.*theta_b(i) theta_b_from_epsilon_zz(i)],'-k','LineWidth',1.1); end
% set(gca,'xtick',1:1:21,'xticklabel',station_names); axe1.YAxis.Exponent = 0;
% legend('from BedMachine','from station-pair $\dot{\epsilon}_{zz,bg}$','Interpreter','latex','Location','SouthEast');
% %for i=1:1:22; text(i+0.2,epsilon_dot_zz_bg(i,1),sprintf('%5.4f',epsilon_dot_zz_bg(i,1))); end
% grid on; xlim([0 22]); %ylim([-0.01 0.01])
% title('bed slope \theta_{b} during 2023/150-160  [ yr^{-1} ]')
% ylabel('$\theta_{b}$  [ degrees ] (positive = bed slopes up against flow)','Interpreter','latex','FontSize',14)
% 
% % print figure
% % cd catalogue_c_dot/30days_48hrs/
% % print(gcf,'-dpng','-r500',sprintf('thetab_2023_150_155_compare_20240121.png')); 
% % cd ../../
% 
% %% scrap figure for colors
% Fig3 = figure(3); clf; 
% set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',[5 2 19 14]);
% axe1 = axes('Position',[0.075 0.70 0.9 0.25],'Box','on','NextPlot','add','XTickLabels',[]);
% axes(axe1)
% h1=plot([168.85 168.85], [0 100],'-'); c1 = get(h1,'Color'); hold on;
% h2=plot([169.85 168.85], [0 100],'-'); c2 = get(h2,'Color');
% h3=plot([170.85 168.85], [0 100],'-'); c3 = get(h3,'Color');
% h4=plot([171.85 168.85], [0 100],'-'); c4 = get(h4,'Color');
% h5=plot([172.85 168.85], [0 100],'-'); c5 = get(h5,'Color');
% h6=plot([173.85 168.85], [0 100],'-'); c6 = get(h6,'Color');
% h7=plot([174.85 168.85], [0 100],'-'); c7 = get(h7,'Color');
% h8=plot([175.85 168.85], [0 100],'-'); c8 = get(h8,'Color');
% h9=plot([176.85 168.85], [0 100],'-'); c9 = get(h9,'Color');
% 
% %% Figures: c_dot_delta_t components and TRACK vertical positions
% % these are for station-paired $\dot{\epsilon}_{zz} H$ 
% for i=1:1:num_sta
% t_plot = [145 175];
% Fig2 = figure(i+10); clf; 
% set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',1.2.*[5 2 19 18]);
% axe1 = axes('Position',[0.08 0.70 0.9 0.27],'Box','on','NextPlot','add','XTickLabels',[]);
% axe2 = axes('Position',[0.08 0.38 0.9 0.27],'Box','on','NextPlot','add','XTickLabels',[]);
% axe3 = axes('Position',[0.08 0.06 0.9 0.27],'Box','on','NextPlot','add');
% 
% axes(axe1)
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_vert_displ,'.-','Color',h1.Color,'MarkerSize',4); hold on;
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t,'.-','Color',h2.Color,'MarkerSize',4); hold on;
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_advective,'.-','Color',h3.Color,'MarkerSize',4);
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_epsilon_zz,'.-','Color',h4.Color,'MarkerSize',4);
% xlim([t_plot(1) t_plot(2)])
% mytitle = ['vertical velocity  [m day$^{-1}$]  (', station_names{i}, ') '];
% title(mytitle,'Interpreter','latex'); 
% ylabel('vertical velocity [m day$^{-1}$]','Interpreter','latex','FontSize',12)
% legend('$w_{s}$','$\dot{c}$','$\frac{u_{s}}{u_{s,bg}}(u_{b,bg}\tan{\theta_{b}})$','$\dot{\epsilon}_{zz} H$',...
%     'Interpreter','latex','Location','NorthWest','NumColumns',2)
% set(gca,'FontName','Avenir','tickdir','in'); grid on
% 
% axes(axe2)
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_vert_displ_cumulative,'.-','Color',h1.Color,'MarkerSize',4); hold on;
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_cumulative,'.-','Color',h2.Color,'MarkerSize',4); hold on;
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_advective_cumulative,'.-','Color',h3.Color,'MarkerSize',4);
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_epsilon_zz_cumulative,'.-','Color',h4.Color,'MarkerSize',4);
% xlim([t_plot(1) t_plot(2)])
% mytitle = ['vertical displacement from 2023/150.5 onwards [m] (', station_names{i}, ') '];
% title(mytitle,'Interpreter','latex'); 
% ylabel('vertical displacement  [m]','Interpreter','latex','FontSize',12)
% legend('$z_{s}$','${c}$','$\frac{u_{s}}{u_{s,bg}}(u_{b,bg}\tan{\theta_{b}})$','$\dot{\epsilon}_{zz} H$',...
%     'Interpreter','latex','Location','NorthWest','NumColumns',2)
% set(gca,'FontName','Avenir','tickdir','in'); grid on
% 
% axes(axe3)
% plot(B2023(i).neu_iu(:,1),B2023(i).neu_iu(:,4),'.','Color',h1.Color,'MarkerSize',2)
% xlim([t_plot(1) t_plot(2)])
% mytitle = ['station vertical position, $z_{s}$   [m above ref. ellipsoid] (', station_names{i}, ') '];
% title(mytitle,'Interpreter','latex'); 
% ylabel('vertical position  [m (relative)]','Interpreter','latex','FontSize',12)
% legend('$z_{s}$, relative vertical position','Interpreter','latex','Location','NorthWest')
% xlabel('DOY in 2023 [ UTC ]'); 
% set(gca,'FontName','Avenir','tickdir','in'); grid on
% 
% % print figure
% print(gcf,'-dpng','-r300',sprintf('catalogue_c_dot/cdot_zz_pairs_2023_150_180_tight_s0_w24_t20_240722_%s.png',station_names{i})); 
% close(Fig2)
% 
% end
% %%
% for i=1:1:num_sta
% t_plot = [165 195];
% Fig2 = figure(i+10); clf; 
% set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',1.2.*[5 2 19 18]);
% axe1 = axes('Position',[0.08 0.70 0.9 0.27],'Box','on','NextPlot','add','XTickLabels',[]);
% axe2 = axes('Position',[0.08 0.38 0.9 0.27],'Box','on','NextPlot','add','XTickLabels',[]);
% axe3 = axes('Position',[0.08 0.06 0.9 0.27],'Box','on','NextPlot','add');
% 
% axes(axe1)
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_vert_displ,'.-','Color',h1.Color,'MarkerSize',4); hold on;
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t,'.-','Color',h2.Color,'MarkerSize',4); hold on;
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_advective,'.-','Color',h3.Color,'MarkerSize',4);
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_epsilon_zz,'.-','Color',h4.Color,'MarkerSize',4);
% xlim([t_plot(1) t_plot(2)])
% mytitle = ['vertical velocity  [m day$^{-1}$]  (', station_names{i}, ') '];
% title(mytitle,'Interpreter','latex'); 
% ylabel('vertical velocity [m day$^{-1}$]','Interpreter','latex','FontSize',12)
% legend('$w_{s}$','$\dot{c}$','$\frac{u_{s}}{u_{s,bg}}(u_{b,bg}\tan{\theta_{b}})$','$\dot{\epsilon}_{zz} H$',...
%     'Interpreter','latex','Location','NorthWest','NumColumns',2)
% set(gca,'FontName','Avenir','tickdir','in'); grid on
% 
% axes(axe2)
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_vert_displ_cumulative,'.-','Color',h1.Color,'MarkerSize',4); hold on;
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_cumulative,'.-','Color',h2.Color,'MarkerSize',4); hold on;
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_advective_cumulative,'.-','Color',h3.Color,'MarkerSize',4);
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_epsilon_zz_cumulative,'.-','Color',h4.Color,'MarkerSize',4);
% xlim([t_plot(1) t_plot(2)])
% mytitle = ['vertical displacement from 2023/150.5 onwards [m] (', station_names{i}, ') '];
% title(mytitle,'Interpreter','latex'); 
% ylabel('vertical displacement  [m]','Interpreter','latex','FontSize',12)
% legend('$z_{s}$','${c}$','$\frac{u_{s}}{u_{s,bg}}(u_{b,bg}\tan{\theta_{b}})$','$\dot{\epsilon}_{zz} H$',...
%     'Interpreter','latex','Location','NorthWest','NumColumns',2)
% set(gca,'FontName','Avenir','tickdir','in'); grid on
% 
% axes(axe3)
% plot(B2023(i).neu_iu(:,1),B2023(i).neu_iu(:,4),'.','Color',h1.Color,'MarkerSize',2)
% xlim([t_plot(1) t_plot(2)])
% mytitle = ['station vertical position, $z_{s}$   [m above ref. ellipsoid] (', station_names{i}, ') '];
% title(mytitle,'Interpreter','latex'); 
% ylabel('vertical position  [m (relative)]','Interpreter','latex','FontSize',12)
% legend('$z_{s}$, relative vertical position','Interpreter','latex','Location','NorthWest')
% xlabel('DOY in 2023 [ UTC ]'); 
% set(gca,'FontName','Avenir','tickdir','in'); grid on
% 
% % print figure
% print(gcf,'-dpng','-r300',sprintf('catalogue_c_dot/cdot_zz_pairs_2023_170_200_tight_s0_w24_t20_240722_%s.png',station_names{i})); 
% close(Fig2)
% 
% end
% 
% for i=1:1:num_sta
% t_plot = [185 215];
% Fig2 = figure(i+10); clf; 
% set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',1.2.*[5 2 19 18]);
% axe1 = axes('Position',[0.08 0.70 0.9 0.27],'Box','on','NextPlot','add','XTickLabels',[]);
% axe2 = axes('Position',[0.08 0.38 0.9 0.27],'Box','on','NextPlot','add','XTickLabels',[]);
% axe3 = axes('Position',[0.08 0.06 0.9 0.27],'Box','on','NextPlot','add');
% 
% axes(axe1)
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_vert_displ,'.-','Color',h1.Color,'MarkerSize',4); hold on;
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t,'.-','Color',h2.Color,'MarkerSize',4); hold on;
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_advective,'.-','Color',h3.Color,'MarkerSize',4);
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_epsilon_zz,'.-','Color',h4.Color,'MarkerSize',4);
% xlim([t_plot(1) t_plot(2)])
% mytitle = ['vertical velocity  [m day$^{-1}$]  (', station_names{i}, ') '];
% title(mytitle,'Interpreter','latex'); 
% ylabel('vertical velocity [m day$^{-1}$]','Interpreter','latex','FontSize',12)
% legend('$w_{s}$','$\dot{c}$','$\frac{u_{s}}{u_{s,bg}}(u_{b,bg}\tan{\theta_{b}})$','$\dot{\epsilon}_{zz} H$',...
%     'Interpreter','latex','Location','NorthWest','NumColumns',2)
% set(gca,'FontName','Avenir','tickdir','in'); grid on
% 
% axes(axe2)
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_vert_displ_cumulative,'.-','Color',h1.Color,'MarkerSize',4); hold on;
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_cumulative,'.-','Color',h2.Color,'MarkerSize',4); hold on;
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_advective_cumulative,'.-','Color',h3.Color,'MarkerSize',4);
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_epsilon_zz_cumulative,'.-','Color',h4.Color,'MarkerSize',4);
% xlim([t_plot(1) t_plot(2)])
% mytitle = ['vertical displacement from 2023/150.5 onwards [m] (', station_names{i}, ') '];
% title(mytitle,'Interpreter','latex'); 
% ylabel('vertical displacement  [m]','Interpreter','latex','FontSize',12)
% legend('$z_{s}$','${c}$','$\frac{u_{s}}{u_{s,bg}}(u_{b,bg}\tan{\theta_{b}})$','$\dot{\epsilon}_{zz} H$',...
%     'Interpreter','latex','Location','NorthWest','NumColumns',2)
% set(gca,'FontName','Avenir','tickdir','in'); grid on
% 
% axes(axe3)
% plot(B2023(i).neu_iu(:,1),B2023(i).neu_iu(:,4),'.','Color',h1.Color,'MarkerSize',2)
% xlim([t_plot(1) t_plot(2)])
% mytitle = ['station vertical position, $z_{s}$   [m above ref. ellipsoid] (', station_names{i}, ') '];
% title(mytitle,'Interpreter','latex'); 
% ylabel('vertical position  [m (relative)]','Interpreter','latex','FontSize',12)
% legend('$z_{s}$, relative vertical position','Interpreter','latex','Location','NorthWest')
% xlabel('DOY in 2023 [ UTC ]'); 
% set(gca,'FontName','Avenir','tickdir','in'); grid on
% 
% % print figure
% print(gcf,'-dpng','-r300',sprintf('catalogue_c_dot/cdot_zz_pairs_2023_190_220_tight_s0_w24_t20_240722_%s.png',station_names{i})); 
% close(Fig2)
% 
% end
% 
% for i=1:1:num_sta
% t_plot = [205 235];
% Fig2 = figure(i+10); clf; 
% set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',1.2.*[5 2 19 18]);
% axe1 = axes('Position',[0.08 0.70 0.9 0.27],'Box','on','NextPlot','add','XTickLabels',[]);
% axe2 = axes('Position',[0.08 0.38 0.9 0.27],'Box','on','NextPlot','add','XTickLabels',[]);
% axe3 = axes('Position',[0.08 0.06 0.9 0.27],'Box','on','NextPlot','add');
% 
% axes(axe1)
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_vert_displ,'.-','Color',h1.Color,'MarkerSize',4); hold on;
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t,'.-','Color',h2.Color,'MarkerSize',4); hold on;
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_advective,'.-','Color',h3.Color,'MarkerSize',4);
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_epsilon_zz,'.-','Color',h4.Color,'MarkerSize',4);
% xlim([t_plot(1) t_plot(2)])
% mytitle = ['vertical velocity  [m day$^{-1}$]  (', station_names{i}, ') '];
% title(mytitle,'Interpreter','latex'); 
% ylabel('vertical velocity [m day$^{-1}$]','Interpreter','latex','FontSize',12)
% legend('$w_{s}$','$\dot{c}$','$\frac{u_{s}}{u_{s,bg}}(u_{b,bg}\tan{\theta_{b}})$','$\dot{\epsilon}_{zz} H$',...
%     'Interpreter','latex','Location','NorthEast','NumColumns',2)
% set(gca,'FontName','Avenir','tickdir','in'); grid on
% 
% axes(axe2)
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_vert_displ_cumulative,'.-','Color',h1.Color,'MarkerSize',4); hold on;
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_cumulative,'.-','Color',h2.Color,'MarkerSize',4); hold on;
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_advective_cumulative,'.-','Color',h3.Color,'MarkerSize',4);
% plot(daily(i).t22, daily_epsilon_zz(i).c_dot_delta_t_epsilon_zz_cumulative,'.-','Color',h4.Color,'MarkerSize',4);
% xlim([t_plot(1) t_plot(2)])
% mytitle = ['vertical displacement from 2023/150.5 onwards [m] (', station_names{i}, ') '];
% title(mytitle,'Interpreter','latex'); 
% ylabel('vertical displacement  [m]','Interpreter','latex','FontSize',12)
% legend('$z_{s}$','${c}$','$\frac{u_{s}}{u_{s,bg}}(u_{b,bg}\tan{\theta_{b}})$','$\dot{\epsilon}_{zz} H$',...
%     'Interpreter','latex','Location','NorthEast','NumColumns',2)
% set(gca,'FontName','Avenir','tickdir','in'); grid on
% 
% axes(axe3)
% plot(B2023(i).neu_iu(:,1),B2023(i).neu_iu(:,4),'.','Color',h1.Color,'MarkerSize',2)
% xlim([t_plot(1) t_plot(2)])
% mytitle = ['station vertical position, $z_{s}$   [m above ref. ellipsoid] (', station_names{i}, ') '];
% title(mytitle,'Interpreter','latex'); 
% ylabel('vertical position  [m (relative)]','Interpreter','latex','FontSize',12)
% legend('$z_{s}$, relative vertical position','Interpreter','latex','Location','NorthEast')
% xlabel('DOY in 2023 [ UTC ]'); 
% set(gca,'FontName','Avenir','tickdir','in'); grid on
% 
% % print figure
% print(gcf,'-dpng','-r300',sprintf('catalogue_c_dot/cdot_zz_pairs_2023_210_240_tight_s0_w24_t20_240722_%s.png',station_names{i})); 
% close(Fig2)
% 
% end
