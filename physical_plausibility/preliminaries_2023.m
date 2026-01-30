% physics-informed cluster chronology 2023 preamble
% 
% load in 2023 preliminary datasets needed for cluster chrono
% 25-04-23 LAS: update to polarstereographic lake locations

%% load north lake geographic files relative to North Lake
origin = [68.72, -49.53];
% polarstereo conversion needed values
radius=6378137.0;    eccen=0.08181919;    lat_true=70;    lon_posy=-45;
addpath('../GNSS_derived/strain_rates_2023/')

%% load prelimiaries for site locations in 2023
% load station names
load('../GNSS_derived/strain_rates_2023/station_names_2023.mat');
num_sta_23 = length(station_names)-1; % number of stations
station_names_short = upper(vertcat(station_names(1,1:18)', station_names(1,20:23)')); % omit SQ33 (#19)
% load north lake geographic files and morlighem bed relative to North Lake
load('../GNSS_derived/strain_rates_2023/polarstereo_stations_2023_short.mat')
% load lake boundaries from FASTER for 2023
%  boundaries and correct dates:
load('../GNSS_derived/strain_rates_2023/environs_lakes_2023B_S1S2MO_250416_noboundaries.mat') %  no boundaries and correct dates
environs_lakes_2023B = environs_lakes; 
X_km = environs_lakes_2023B.X_km; 
Y_km = environs_lakes_2023B.Y_km;
Surface_elev = environs_lakes_2023B.S;

% within subglacial ROI
ROI_relevant_xv = [-15 60 104 104 15 -10 -33 -33 -15]; % km in nevis-model space
ROI_relevant_yv = [15 4 -8 -50 -50 -41 -33 -1 15];
for i=1:1:length(environs_lakes.H)
    if inpolygon(environs_lakes.X_km(i),environs_lakes.Y_km(i),...
            ROI_relevant_xv,ROI_relevant_yv) == 0
        % throw a NaN
        environs_lakes.H(i) = NaN;
        environs_lakes.S(i) = NaN;
        environs_lakes.X_km(i) = NaN;
        environs_lakes.Y_km(i) = NaN;
        environs_lakes.lat(i) = NaN;
        environs_lakes.lon(i) = NaN;
        environs_lakes.drainage_type_num(i) = NaN;
        environs_lakes.laketypeing_dates(i,1:7) = [NaN, NaN, NaN, NaN, NaN, NaN, NaN];
    else
    end
end
% redefine lake_typeing to just include lakes within nevis domain
lake_typeing = environs_lakes.laketypeing_dates; % after throwing out lakes on the nevis domain border
type_num = lake_typeing(:,7); % what type of drainage 
lake_num = length(type_num); % number of lakes
environs_lakes_2023 = environs_lakes; 

% polygon ROI check
figure; plot(ROI_relevant_xv,ROI_relevant_yv)
hold on; plot(environs_lakes.X_km,environs_lakes.Y_km,'o')

% load bedmap
load('../GNSS_derived/strain_rates_2022/cmapland2.mat'); % bed topo colormap  
load('../GNSS_derived/strain_rates_2022/BMv5_for_nevis_catchment.mat'); % origin = M1 moulin 

%% load in datasets for nevis model and drainage classifier from S1/2
% drainage type-ing
% 1 = HF; 2 = moulin; 3 = overspill; 4 = crevasses; 5 = frozen. 
% 6= slush; 7 = stream; 9 = throw out/inconclusive.

% load daily_lake_classifier timeseries (all possible S1+S2 images have been used)
% lakes are only around when they are voluminous to HF
load('../GNSS_derived/strain_rates_2023/daily_lake_HFpossible_classifier_S1S2MO_2023_240903.mat')
daily_lake_classifier = daily_lake_HFpossible_classifier;
daily_lake_classifier(44,189)=1; % still-there HF, on inner S2 image bound
daily_lake_classifier(66,189)=1; % still-there HF, on inner S2 image bound
daily_lake_classifier(103,189)=1; % still-there HF, on inner S2 image bound
daily_lake_classifier(131,189)=1; % still-there HF, on inner S2 image bound
% within subglacial ROI:
for i=1:1:length(environs_lakes_2023.S)
    if isnan(environs_lakes_2023.S(i))
        % throw a NaN
        daily_lake_classifier(i,:) = NaN;
    else
    end
end

% lake locations
load('../GNSS_derived/strain_rates_2023/max_extent_indices_2023_GEOD.mat');

% bedmap in nevis formation
load('../GNSS_derived/strain_rates_2022/BMv5_for_nevis_catchment_noSK.mat'); % bed machine
H = BMv5_for_nevis_catchment_noSK.S_m - BMv5_for_nevis_catchment_noSK.B_m; % ice thickness
S = BMv5_for_nevis_catchment_noSK.S_m; % ice surface elevation

% load nevis discharge q outputs
    addpath('../nevis_lakes_cluster/nevis/')
    addpath('../nevis_lakes_cluster/nevis_lakesix_2023_noSK_300m_ub/')
    % load initial timestep
    fn = '../nevis_lakes_cluster/nevis_lakesix_2023_noSK_300m_ub';
    if exist([fn,'/0000']), load([fn,'/0000']); else load(fn); end
    load([fn,'/0001']);
    if isfield(oo,'reversey'), reversey = oo.reversey; else reversey = 0; end
    
    % extract variables
    if isfield(vv,'nbdy'), gg = nevis_label(gg,vv.nbdy); aa.phi = aa.phi_a(gg.nbdy); end % re-mask
    aa = nevis_inputs(vv.t,aa,pp,gg,oo);
    oo.evaluate_variables = 1; [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
    vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
    nevis_unpack(aa,gg,vv2);
    
    %get rid of points outside domain
    nx(gg.nout) = NaN;
    ex(gg.eout) = NaN;
    fx(gg.fout) = NaN;
    cx(gg.cout) = NaN;

    % x and y grid:
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny; % km

    % ice thickness
    HH = reshape(ps.z*H,fI,eJ); 
    SS = HH+ps.z*reshape(b,gg.nI,gg.nJ);  
        
    %boundary curve
    if ~isempty(gg.n1)
    x_out = gg.nx(gg.n1); y_out = gg.ny(gg.n1);
    tmp = nevis_orderboundary(x_out,y_out); x_out = x_out(tmp); y_out = y_out(tmp); % reorder to follow boundary
    else x_out = []; y_out = [];
    end

% ID_lakes_nevis_grid
xx_vec = reshape(xx,601*264,1); 
yy_vec = reshape(yy,601*264,1);
for i=1:1:length(X_km)
   % distance of x (and y) point from x (and y) location in RACMO grid
   offset_x =  X_km(1,i) - xx_vec; % km 
   offset_y =  Y_km(1,i) - yy_vec; % km
   % distance of XY point from all XY points in the nevis grid
   offset_xy = sqrt((offset_x.^2)+(offset_y.^2));
   [dist, ID_xy] = min(offset_xy); % find minimum offset: index and value
   ID_lakes_nevis_grid(i,1:2) = [ID_xy, dist]; % keep index and kms away from RACMO grid location
end

% load in nevis outputs for Q_out, basal melt, and E across domain
load('../nevis_lakes_cluster/nevis_lakesix_2023_noSK_300m_ub_unique.mat')

%% winter stresses from MEaSUREs Annual or Monthly compilation Sentinel-1 SAR
load('winter-stress-datasets/stresses_nlake_2023Annual_S12_250219.mat');
stresses_XB = stresses_nlake_2023Annual_S12.XB_strain_grid_sub; % [ m ]
stresses_YB = stresses_nlake_2023Annual_S12.YB_strain_grid_sub; % [ m ]
stresses_sig1_winter = stresses_nlake_2023Annual_S12.Smax; % [ Pa ]
flow_angle = stresses_nlake_2023Annual_S12.flow_angle; % [ degrees ]

% calculate sig1_winter at lake locations:
stresses_XB_vec = reshape(stresses_XB(:,1:end-1), 201*225, 1);
stresses_YB_vec = reshape(stresses_YB(:,1:end-1), 201*225, 1);
stresses_sig1_winter_vec = reshape(stresses_sig1_winter, 201*225, 1);
for i=1:1:length(X_km)
   % distance of x (and y) point from x (and y) location in RACMO grid
   offset_x =  X_km(1,i) - stresses_XB_vec./1e3; % km 
   offset_y =  Y_km(1,i) - stresses_YB_vec./1e3; % km
   % distance of XY point from all XY points in the RACMO grid
   offset_xy = sqrt((offset_x.^2)+(offset_y.^2));
   [dist, ID_xy] = min(offset_xy); % find minimum offset: index and value
   ID_stresses_sig1_winter(i,1:2) = [ID_xy, dist]; % keep index and kms away from RACMO grid location
end