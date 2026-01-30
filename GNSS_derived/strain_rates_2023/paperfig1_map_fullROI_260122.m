% strain choices MAP figure
% a figure that shows choices made for calculating c_dot at GNSS sites both
% in 2022 and 2023, also includes FASTER lake outlines :D

% LAS 24-08-19: first go
% LAS 24-08-21: colors VIVID
% LAS 24-08-23: moulin drainage times back in, full ROI shown
% LAS 24-08-30: S1+S2 drainage dates :)
% 24-09-03 LAS: time of drainage plotted as DOY values of 24-hr periods, denote freezers
% 24-09-06 LAS: ROI polygon inspired by subglacial hydro catchments
% 25-04-16 LAS: consistent XY-grid for FASTER lake boundaries
% 25-06-04 LAS: ice-flow direction in map-view panels
% 25-09-26 LAS: add in Natalie Turner's moulins (2022)
% 25-10-18 LAS: in-ROI statistics
% 26-01-22 LAS: +longterm ELA 

close all; clear all

%% LOCATION MAP %% THIS BEDMACHINE SUBSET IS TOO LARGE FOR FIGSHARE
% load bedmap
% load BM3_westcoast_inland.mat % westcoast bedmachine 3
% mask=double(BM3_westcoast.Mask);
% %0 = ocean, 1 = ice-free land, 2 = grounded ice, 3 = floating ice
% mask(mask==0)=NaN; mask(mask==1)=NaN; mask(mask==3)=NaN; mask(mask==2)=1;
% Surface_Ice = mask.*(double(BM3_westcoast.Thickness)+double(BM3_westcoast.Bed));
% Surface = double(BM3_westcoast.Surface);
% 
% BM3_lon = BM3_westcoast.lon; BM3_lat = BM3_westcoast.lat;
% 
% mask2=double(BM3_westcoast.Mask);
% mask2(mask2==2)=NaN; mask2(mask2==3)=NaN; mask2(mask2==0)=1;
% Surface_Not_Ice = mask2.*(double(BM3_westcoast.Surface));
% Bed_Not_Ice = mask2.*(double(BM3_westcoast.Bed));
% 
% lat_vec = reshape(BM3_lat, 2334*4335, 1);
% lon_vec = reshape(BM3_lon, 2334*4335, 1);
% Surface_Not_Ice_vec = reshape(Surface_Not_Ice, 2334*4335, 1);
% Bed_vec = reshape(BM3_westcoast.Bed, 2334*4335, 1);
% Bed_Not_Ice_vec = reshape(Bed_Not_Ice, 2334*4335, 1);
% Surface_Ice_vec = reshape(Surface_Ice, 2334*4335, 1);
% 
% surface_nan = horzcat(lon_vec, lat_vec, Surface_Ice_vec);
% surface_nan(any(isnan(surface_nan), 2), :) = [];
% clear BM3_westcoast

%% load prelimiaries for site locations 2022
% XY relative to moulin M1 in Stevens et al. (2015)
origin = [68.72, -49.53]; % M1 moulin
% polarstereo conversion needed values
radius=6378137.0;    eccen=0.08181919;    lat_true=70;    lon_posy=-45;
[moulin_x,moulin_y] = polarstereo_fwd(origin(1),origin(2),radius,eccen,lat_true,lon_posy); % [ m ]
moulin_x_km = moulin_x./1e3; moulin_y_km = moulin_y./1e3; % UTM in KM of moulin
% load station names
station_names_2022 = load('../strain_rates_2022/station_names.mat');
num_sta_22 = length(station_names_2022.station_names)-2; % number of stations
station_names_2022_short = upper(vertcat(station_names_2022.station_names(1,1:17)',...
    station_names_2022.station_names(1,20:22)')); % omit SQ32, SQ33 (#18, #19)
% load site locations 2022
load('../strain_rates_2022/apcoords_lle_2022.mat');
    lats=apcoords_lle_2022(1,:); lons=apcoords_lle_2022(2,:)-360; hs=apcoords_lle_2022(3,:); llh=[lats; lons; hs];
    [xy_sta_22(:,1), xy_sta_22(:,2)] = polarstereo_fwd(lats,lons,radius,eccen,lat_true,lon_posy); % [ m ]
    xy_sta_22(:,1) = (xy_sta_22(:,1)-moulin_x)./1e3; % [ km ] 
    xy_sta_22(:,2) = (xy_sta_22(:,2)-moulin_y)./1e3; % [ km ] 
    xy_sta_22_short = vertcat(xy_sta_22(1:17,:), xy_sta_22(19:end,:)); % omit SQ33 (#18)
save polarstereo_stations_2022_short.mat xy_sta_22_short
% load lake boundaries, centre points, drainage dates from FASTER for 2022
load('../strain_rates_2022/environs_lakes_2022B_250416.mat') % drainage dates & boundaries
environs_lakes_2022 = environs_lakes; 
% load bedmap
load('../strain_rates_2022/BMv5_for_nevis_catchment.mat'); % origin = M1 moulin 
load('../strain_rates_2022/cmapland2.mat'); % bed topo colormap  
% load wintertime velocity map
load('../../physical_plausibility/winter-stress-datasets/stresses_nlake_2022Annual_S12_250219.mat')
% inerpret to nevis grid for nevis ice mask
flow_angle_interp_2022 = interp2(stresses_nlake_2022Annual_S12.XB_strain_grid_sub,...
    stresses_nlake_2022Annual_S12.YB_strain_grid_sub,stresses_nlake_2022Annual_S12.flow_angle,...
    BMv5_for_nevis_catchment.X_km,BMv5_for_nevis_catchment.Y_km);
mask = isfinite(BMv5_for_nevis_catchment.B_km);
mask = double(mask); mask(mask==0)=NaN;
flow_angle_interp_2022 = mask.*flow_angle_interp_2022;
% load Natalie Turner's moulins
NT = shaperead('../strain_rates_2022/NatalieTurner_Moulins_2022/moulins.shp');
[moulins_x,moulins_y] = polarstereo_fwd([NT.Y],[NT.X],radius,eccen,lat_true,lon_posy); % [ m ]
moulins_x_2022 = (moulins_x-moulin_x)./1e3; moulins_y_2022 = (moulins_y-moulin_y)./1e3;
% load Natalie Turner's connections to lakes and moulins
load('../strain_rates_2022/NatalieTurner_Moulins_2022/laketypeing_2022_250930_result.mat');

%% load prelimiaries for site locations 2023
% load station names
load station_names_2023.mat
num_sta = length(station_names)-1; % number of stations
station_names_short = upper(vertcat(station_names(1,1:18)', station_names(1,20:23)')); % omit SQ33 (#19)
% load site locations 2023
load apcoords_lle_2023.mat
    lats=apcoords_lle_2023(1,:); lons=apcoords_lle_2023(2,:); 
    [xy_sta_23(:,1), xy_sta_23(:,2)] = polarstereo_fwd(lats,lons,radius,eccen,lat_true,lon_posy); % [ m ]
    xy_sta_23 = vertcat(xy_sta_23(22:23,:), xy_sta_23(1:21,:)); % [ m ] 
    xy_sta_23(:,1) = (xy_sta_23(:,1)-moulin_x)./1e3; % [ km ] 
    xy_sta_23(:,2) = (xy_sta_23(:,2)-moulin_y)./1e3; % [ km ] 
    xy_sta_23_short = vertcat(xy_sta_23(1:18,:), xy_sta_23(20:23,:)); % omit SQ33 (#19)
save polarstereo_stations_2023_short.mat xy_sta_23_short
% load lake centre points, boundaries from FASTER for 2023
%  boundaries and correct dates:
load('environs_lakes_2023B_S1S2MO_251003.mat') %  boundaries and correct dates
environs_lakes_2023B = environs_lakes; 
% load wintertime velocity map
load('../../physical_plausibility/winter-stress-datasets/stresses_nlake_2023Annual_S12_250219.mat')
% inerpret to nevis grid for nevis ice mask
flow_angle_interp_2023 = interp2(stresses_nlake_2023Annual_S12.XB_strain_grid_sub,...
    stresses_nlake_2023Annual_S12.YB_strain_grid_sub,stresses_nlake_2023Annual_S12.flow_angle,...
    BMv5_for_nevis_catchment.X_km,BMv5_for_nevis_catchment.Y_km);
flow_angle_interp_2023 = mask.*flow_angle_interp_2023;
% load Laura's 2023 moulins
NT_23 = shaperead('../strain_rates_2022/NatalieTurner_Moulins_2022/moulins_2023.shp');
[moulins23_x,moulins23_y] = polarstereo_fwd([NT_23.Y],[NT_23.X],radius,eccen,lat_true,lon_posy); % [ m ]
moulins_x_2023 = (moulins23_x-moulin_x)./1e3; moulins_y_2023 = (moulins23_y-moulin_y)./1e3;
% load Laura's connections to lakes and moulins
load('../strain_rates_2022/NatalieTurner_Moulins_2022/laketypeing_2023_251002_result.mat');

%% RACMO longterm ELA 1958–2019 Brice Noël et a. (2019) Sci Advances
ELA = shaperead('../strain_rates_2022/NatalieTurner_Moulins_2022/ELA-longterm/RACMO_1km_SMBZERO_contour.shp');
% last entry is the ELA along CWS to S GrIS: 1074
[ELA(1074).lat, ELA(1074).lon] = polarstereo_inv([ELA(1074).X]', [ELA(1074).Y]',radius,eccen,lat_true,lon_posy);
ELA(1074).X_plot = (([ELA(1074).X]')-moulin_x)./1e3;
ELA(1074).Y_plot = (([ELA(1074).Y]')-moulin_y)./1e3;

%% location map + 2022 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf; 
set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',1.34.*[1 0 18.3 22]);
axe0 = axes('Position',[0.05 0.825 0.35 0.16],'Box','on','NextPlot','add');
axe1 = axes('Position',[0.05 0.425 0.945 0.40],'Box','on','NextPlot','add','XTickLabel',[]);
axe2 = axes('Position',[0.05 0.03 0.945 0.40],'Box','on','NextPlot','add');

% plotting preliminaries
fontsize = 9; 
TriangleSize = 9;
DiamondSize = 0.5;
linewidth_strain = 0.9;
mlat = 68.8;
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
moulin_green = [0.4 0.8 0.8];
lake_colors = vertcat(goldenrod, moulin_green, lake_blue(1:3), icy_plum);
% full ROI
plot_x = [-34 104];
plot_y = [-51.5 16.5];

% location map
axes(axe0)
text(-52.70, 69.35,'a','FontWeight','bold','FontSize',fontsize+3); 
hold on; hold all;

% % BedMap plotting:
% sss=scatter(lon_vec,lat_vec,2,Bed_Not_Ice_vec,'s','fill');
% set(sss,'MarkerFaceAlpha',0.8); 
% load oceanbone.mat
% demcmap([-1000 1000]); colormap(oceanbone);
% hold all
% vv=200:200:1800;
% [C,Ch]=contour(BM3_lon, BM3_lat, Surface_Ice,vv);
% set(Ch,'LineColor',[0.5 0.5 0.5],'LineWidth',1.05);
% clabel(C,Ch,'FontSize',fontsize-1,'LabelSpacing',175,'Color',[0.5 0.5 0.5])

% lakes: FASTER S2 2022
[HF_lakes_to_plot_22] = find(environs_lakes_2022.laketypeing_dates(:,7)==1); % just HF lakes   
[MU_lakes_to_plot_22] = find(environs_lakes_2022.laketypeing_dates(:,7)==2); % just moulin-draining lakes 
[OS_lakes_to_plot_22] = find(environs_lakes_2022.laketypeing_dates(:,7)==3); % just overspill-draining lakes 
[NE_lakes_to_plot_22] = find(environs_lakes_2022.laketypeing_dates(:,7)==5); % just no-exit freezers
[lakes_to_plot_22] = vertcat(HF_lakes_to_plot_22, MU_lakes_to_plot_22, OS_lakes_to_plot_22, NE_lakes_to_plot_22); % all good lakes   
% plotting in order for legend
fill3(environs_lakes_2022.boundaries(lakes_to_plot_22(1)).lon(:,1), ...
        environs_lakes_2022.boundaries(lakes_to_plot_22(1)).lat(:,1), ...
        3000.*ones(length(environs_lakes_2022.boundaries(lakes_to_plot_22(1)).lon(:,1)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',lake_blue(1:3));
fill3(environs_lakes_2022.boundaries(HF_lakes_to_plot_22(1)).lon(:,1), ...
        environs_lakes_2022.boundaries(HF_lakes_to_plot_22(1)).lat(:,1), ...
        3000.*ones(length(environs_lakes_2022.boundaries(HF_lakes_to_plot_22(1)).lon(:,1)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',goldenrod);
fill3(environs_lakes_2022.boundaries(MU_lakes_to_plot_22(1)).lon(:,1), ...
        environs_lakes_2022.boundaries(MU_lakes_to_plot_22(1)).lat(:,1), ...
        3000.*ones(length(environs_lakes_2022.boundaries(MU_lakes_to_plot_22(1)).lon(:,1)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',[0.4 0.8 0.8]);
fill3(environs_lakes_2022.boundaries(NE_lakes_to_plot_22(1)).lon(:,1), ...
        environs_lakes_2022.boundaries(NE_lakes_to_plot_22(1)).lat(:,1), ...
        3000.*ones(length(environs_lakes_2022.boundaries(NE_lakes_to_plot_22(1)).lon(:,1)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',plum(1:3));
% fill those lakes, girl
for i=1:1:length(lakes_to_plot_22)
    fill3(environs_lakes_2022.boundaries(lakes_to_plot_22(i)).lon(:,1), ...
        environs_lakes_2022.boundaries(lakes_to_plot_22(i)).lat(:,1), ...
        3000.*ones(length(environs_lakes_2022.boundaries(lakes_to_plot_22(i)).lon(:,1)),1),... 
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',lake_blue(1:3));
end
for i=1:1:length(HF_lakes_to_plot_22)
    fill3(environs_lakes_2022.boundaries(HF_lakes_to_plot_22(i)).lon(:,1), ...
        environs_lakes_2022.boundaries(HF_lakes_to_plot_22(i)).lat(:,1), ...
        3000.*ones(length(environs_lakes_2022.boundaries(HF_lakes_to_plot_22(i)).lon(:,1)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',goldenrod);
end
for i=1:1:length(MU_lakes_to_plot_22)
    fill3(environs_lakes_2022.boundaries(MU_lakes_to_plot_22(i)).lon(:,1), ...
        environs_lakes_2022.boundaries(MU_lakes_to_plot_22(i)).lat(:,1), ...
         3000.*ones(length(environs_lakes_2022.boundaries(MU_lakes_to_plot_22(i)).lon(:,1)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',[0.4 0.8 0.8]);
end
for i=1:1:length(NE_lakes_to_plot_22)
    fill3(environs_lakes_2022.boundaries(NE_lakes_to_plot_22(i)).lon(:,1), ...
        environs_lakes_2022.boundaries(NE_lakes_to_plot_22(i)).lat(:,1), ...
        3000.*ones(length(environs_lakes_2022.boundaries(NE_lakes_to_plot_22(i)).lon(:,1)),1),...
        'o','MarkerSize',1,'EdgeColor','none','FaceColor',plum(1:3));
end

% ROI that's relevant to the GNSS array; after FASTER locations update (2025-April):
ROI_relevant_xv = [-15 60 104 104 15 -10 -33 -33 -15]; % km in nevis-model space
ROI_relevant_yv = [15 4 -8 -50 -50 -41 -33 -1 15];
ROI_relevant_250417 = vertcat(ROI_relevant_xv, ROI_relevant_yv);
save ROI_relevant_250417.mat ROI_relevant_250417

% longterm ELA 
plot(ELA(1074).lon, ELA(1074).lat,'-','LineWidth',1.5,'Color',[0.2 0.6 0.4]); % nice green

% polarstereo conversion needed values
    radius=6378137.0;    eccen=0.08181919;    lat_true=70;    lon_posy=-45;  
    origin = [68.72, -49.53]; % North Lake M1 moulin in 2011
[moulin_x,moulin_y] = polarstereo_fwd(origin(1),origin(2),radius,eccen,lat_true,lon_posy); % [ m ]
moulin_x_km = moulin_x./1e3; moulin_y_km = moulin_y./1e3; % UTM in KM of moulin
[ROI_relevant_lat, ROI_relevant_lon] = polarstereo_inv((ROI_relevant_xv+moulin_x_km)*1e3, ...
    (ROI_relevant_yv+moulin_y_km)*1e3,radius, eccen, lat_true, lon_posy);
plot(ROI_relevant_lon,ROI_relevant_lat,'-','LineWidth',1.2,'Color',[0.6 0.2 0.0]);
text(-47.95,68.83,10e3,'ROI','FontName','AvenirBold','FontAngle','Italic',...
    'FontSize',fontsize+1,'Rotation',-14,'Color',[0.6 0.2 0.0])
text(-49.4, 69.20, 10e3, 'Sermeq Kujalleq','FontName','AvenirBold','FontAngle','Italic',...
    'BackgroundColor',[0.9 0.9 0.9],...
    'FontSize',fontsize,'Rotation',0,'Color',lake_blue(1:3))
text(-47.15, 69.0, 10e3, 'ELA','FontName','AvenirBold','FontAngle','Italic',...
    'FontSize',fontsize+1,'Rotation',-75,'Color',[0.2 0.6 0.4]);

% scalebar
plot([-51.4 -50.785]+4.2,[68.65 68.65]+0.50,'k-','LineWidth',3);
text(-51.3+4.05, 68.7+0.52, '25 km','Color','k','FontWeight','bold','FontSize',fontsize+1);
% year of lakes catalogue
text(-51.95,69.23,'2022','FontWeight','bold','FontSize',fontsize+2,...
    'Color','w')
xlim([-52 -46.4]); ylim([68.20 69.30]);
set(gca,'dataaspectratio',[1,cos(mlat*pi/180),1],'LineWidth',0.8,...
    'ticklength',[0.015 0.015], 'tickdir','in','FontName','Avenir','Layer','top',...
    'xaxislocation','top','YTick',68.5:0.5:69,'YTickLabel',{'68.5 ˚N'; '69.0 ˚N'},...
    'XTick',-51:1:-47,'FontSize',fontsize,'XTickLabel',[]); % 'XTickLabel',{'51 ˚W';'50 ˚W';'49 ˚W';'48 ˚W';'47 ˚W'});
axe0.ClippingStyle = "rectangle";  
ax=gca; ax.Layer = 'top';
box on; 
text((-51:1:-47)-0.25,[69.40,69.40,69.40,69.40,69.40]-0.04,...
    {'51 ˚W';'50 ˚W';'49 ˚W';'48 ˚W';'47 ˚W'},'FontName','Avenir','FontSize',fontsize)

% mini Green Greenland
axe13 = axes('Position',[0.435 0.84 0.15 0.1833],'NextPlot','add','Box','off',...
    'XTickLabel',[],'YTickLabel',[]);
axes(axe13)
ax = worldmap('greenland');
greenland = shaperead('landareas', 'UseGeoCoords', true,...
  'Selector',{@(name) strcmp(name,'Greenland'), 'Name'});
patchm(greenland.Lat,greenland.Lon,'FaceColor',[0.2 0.6 0.4],'EdgeColor',lake_blue(1:3),'LineWidth',0.8); hold on;
fillm([68.10 69.3 69.3 68.1],[-52 -52 -46.4 -46.4],'w','LineWidth',0.8,'FaceColor',goldenrod,'EdgeColor',goldenrod)
axis off; framem on; gridm on; mlabel off; plabel off;
setm(gca,'FLineWidth',1.1,'FLatLimit',[58 85],'FLonLimit',[-35 35],...
    'FEdgeColor',lake_blue(1:3))

%% 2022 
axes(axe1) 
ylabel('North-South [ km ]','FontSize',fontsize);
text(-40, 14.5, 'b','FontWeight','bold','FontSize',fontsize+3);

hold on;
% ice-sheet bed contour map, colorbar
surf(BMv5_for_nevis_catchment.X_km/1e3,BMv5_for_nevis_catchment.Y_km./1e3,...
   BMv5_for_nevis_catchment.B_km,'EdgeColor','None','Facealpha',0.1); % plot surface below station points
view([0 90])
clim([-500 500]); colormap(cmap)

t2=colorbar('south');
set(t2,'XTick',[-500,-250,0,250,500],'TickDirection','out'); 
hold all
set(get(t2,'xlabel'),'String',{'Bed Elevation  [ m a.s.l. ]'},'FontSize',fontsize,'Fontname','Avenir');
set(t2, 'Position', [0.425 0.83 .17 0.005],'XTick',[-500,-250,0,250,500],...
   'XTickLabel',[-500,-250,0,250,500],'xaxislocation','top');

hold all
% surface ice sheet contours
surf_contours = 0:100:1900; % m a.s.l.
[C, h] = contour(BMv5_for_nevis_catchment.X_km./1e3, BMv5_for_nevis_catchment.Y_km./1e3, ...
    BMv5_for_nevis_catchment.S_km, surf_contours,'Color',ice_sheet_contours,'LineWidth',0.8);
clabel(C,h,'FontSize',7,'Color',ice_sheet_contours,'FontName','Avenir','LabelSpacing',400,'Clipping','on')

% ice-flow direction (for legend)
plot([200 200],[200 200],'-','LineWidth',1.1,'Color',[0.5 0.5 0.5]);
% [0.6 0.2 0.0] % this is a nice, brown colour, but it's too much

% main array GNSS
plot3(xy_sta_22_short(:,1),xy_sta_22_short(:,2),300.*ones(num_sta_22,1),'k^',...
    'MarkerSize',TriangleSize-3,'MarkerFaceColor','w','markerEdgeColor','k',...
    'LineWidth',0.8);

% Natalie Turner's moulins, JUST IN ROI
for i=1:length(moulins_y_2022); IN_M(i) = inpolygon(moulins_x_2022(i),moulins_y_2022(i),ROI_relevant_xv, ROI_relevant_yv); end
moulins_x_2022_INPOLY = moulins_x_2022(IN_M); 
moulins_y_2022_INPOLY = moulins_y_2022(IN_M);
plot(moulins_x_2022_INPOLY(10), moulins_y_2022_INPOLY(10),'ko',...
     'MarkerSize',TriangleSize-7,'MarkerFaceColor','w','markerEdgeColor','k',...
     'LineWidth',0.5);

% lakes: FASTER S2 2022
% [HF_lakes_to_plot_22] = find(environs_lakes_2022.laketypeing_dates(:,7)==1); % just HF lakes   
% [MU_lakes_to_plot_22] = find(environs_lakes_2022.laketypeing_dates(:,7)==2); % just moulin-draining lakes 
% [OS_lakes_to_plot_22] = find(environs_lakes_2022.laketypeing_dates(:,7)==3); % just overspill-draining lakes 
% [lakes_to_plot_22] = find(~isnan(environs_lakes_2022.laketypeing_dates(:,7))); % all good lakes   
% [NE_lakes_to_plot_22] = find(environs_lakes_2022.laketypeing_dates(:,7)==5); % just no-exit freezers

% lakes: FASTER S2 2022 -- IN_ROI_POLYGON
for i=1:length(HF_lakes_to_plot_22)
    IN_HF(i,1) = inpolygon(environs_lakes_2022.X_km(HF_lakes_to_plot_22(i)),environs_lakes_2022.Y_km(HF_lakes_to_plot_22(i)),...
        ROI_relevant_xv, ROI_relevant_yv);
end
[HF_lakes_to_plot_22_INPOLY] = HF_lakes_to_plot_22(IN_HF);
for i=1:length(MU_lakes_to_plot_22)
    IN_MU(i,1) = inpolygon(environs_lakes_2022.X_km(MU_lakes_to_plot_22(i)),environs_lakes_2022.Y_km(MU_lakes_to_plot_22(i)),...
        ROI_relevant_xv, ROI_relevant_yv);
end
[MU_lakes_to_plot_22_INPOLY] = MU_lakes_to_plot_22(IN_MU);
for i=1:length(OS_lakes_to_plot_22)
    IN_OS(i,1) = inpolygon(environs_lakes_2022.X_km(OS_lakes_to_plot_22(i)),environs_lakes_2022.Y_km(OS_lakes_to_plot_22(i)),...
        ROI_relevant_xv, ROI_relevant_yv);
end
[OS_lakes_to_plot_22_INPOLY] = OS_lakes_to_plot_22(IN_OS);
for i=1:length(lakes_to_plot_22)
    IN(i,1) = inpolygon(environs_lakes_2022.X_km(lakes_to_plot_22(i)),environs_lakes_2022.Y_km(lakes_to_plot_22(i)),...
        ROI_relevant_xv, ROI_relevant_yv);
end
[lakes_to_plot_22_INPOLY] = lakes_to_plot_22(IN);
for i=1:length(NE_lakes_to_plot_22)
    NE_IN(i,1) = inpolygon(environs_lakes_2022.X_km(NE_lakes_to_plot_22(i)),environs_lakes_2022.Y_km(NE_lakes_to_plot_22(i)),...
        ROI_relevant_xv, ROI_relevant_yv);
end
[NE_lakes_to_plot_22_INPOLY] = NE_lakes_to_plot_22(NE_IN);

% count up overspills-->other lakes; overspills-->terminal moulins
ROI_area = polyarea(ROI_relevant_xv,ROI_relevant_yv)
in_ROI_lakes_2022 = length(lakes_to_plot_22_INPOLY)
lake_INPOLY_density = length(lakes_to_plot_22_INPOLY)./ROI_area
HF_in_ROI_lakes_2022 = length(HF_lakes_to_plot_22_INPOLY)
MU_in_ROI_lakes_2022 = length(MU_lakes_to_plot_22_INPOLY)
OS_in_ROI_lakes_2022 = length(OS_lakes_to_plot_22_INPOLY)
NE_in_ROI_lakes_2022 = length(NE_lakes_to_plot_22_INPOLY)
OS_into_lake_2022 = sum(strcmp(result{OS_lakes_to_plot_22_INPOLY,3}, 'L'))
OS_into_independent_moulin_2022 = sum(strcmp(result{OS_lakes_to_plot_22_INPOLY,3}, 'M'))
OS_into_NaN_2022 = sum(strcmp(result{OS_lakes_to_plot_22_INPOLY,3}, 'N'))

% plotting in order for legend -- wierdness in ordering due to funky legend
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
for j=21
    i = OS_lakes_to_plot_22_INPOLY(j);
    if strcmp(result{i,3}, 'L')
    plot3([environs_lakes_2022.X_km(i) environs_lakes_2022.X_km(result{i,4})],...
          [environs_lakes_2022.Y_km(i) environs_lakes_2022.Y_km(result{i,4})],[300 300],...
          '-','LineWidth',0.6,'Color',[0 0.55 1])
    else end
end

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
        'MarkerSize',TriangleSize-7,'MarkerFaceColor','w','markerEdgeColor','k',...
        'LineWidth',0.5);
    else end
end

% ice-flow direction
span = 33; tick_scale = 1.0;
quiver3(BMv5_for_nevis_catchment.X_km(1:span:end,1:span:end)./1e3,...
    BMv5_for_nevis_catchment.Y_km(1:span:end,1:span:end)./1e3,...
    1e3.*ones(size(BMv5_for_nevis_catchment.X_km(1:span:end,1:span:end))),... 
    tick_scale.*cosd(flow_angle_interp_2022(1:span:end,1:span:end)),...
    tick_scale.*sind(flow_angle_interp_2022(1:span:end,1:span:end)),...
    0.*sind(flow_angle_interp_2022(1:span:end,1:span:end)),...
    'LineWidth',0.7,'ShowArrowHead','off','Color',[0.5 0.5 0.5],'AutoScale','off');
% [0.6 0.2 0.0] % this is a nice, brown colour, but it's too much

% main array GNSS
plot3(xy_sta_22_short(:,1),xy_sta_22_short(:,2),700.*ones(num_sta_22,1),'k^',...
    'MarkerSize',TriangleSize-3,'MarkerFaceColor','w','markerEdgeColor','k',...
    'LineWidth',0.8);

% ROI that's relevant to the GNSS array
plot(ROI_relevant_xv,ROI_relevant_yv,'-','LineWidth',1.5,'Color',[0.6 0.2 0.0]);
text(81.5,-48,10e3,'Region of Interest (ROI)','FontName','AvenirBold','FontAngle','Italic',...
    'FontSize',fontsize+2,'Color',[0.6 0.2 0.0],'HorizontalAlignment','center')

% longterm ELA 
plot(ELA(1074).X_plot, ELA(1074).Y_plot,...
    '-','LineWidth',1.5,'Color',[0.2 0.6 0.4]); % nice green
text(95.5, -33, 10e3, 'ELA','FontName','AvenirBold','FontAngle','Italic',...
    'FontSize',fontsize+1,'Rotation',-90,'Color',[0.2 0.6 0.4]);

% all lakes
for i=1:1:length(lakes_to_plot_22)
    fill3(environs_lakes_2022.boundaries(lakes_to_plot_22(i)).XY_km_local(:,1), ...
        environs_lakes_2022.boundaries(lakes_to_plot_22(i)).XY_km_local(:,2), ...
        3000.*ones(length(environs_lakes_2022.boundaries(lakes_to_plot_22(i)).XY_km_local(:,2)),1),... 
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',lake_blue(1:3));
end
% HF events
for i=1:1:length(HF_lakes_to_plot_22)
    fill3(environs_lakes_2022.boundaries(HF_lakes_to_plot_22(i)).XY_km_local(:,1), ...
        environs_lakes_2022.boundaries(HF_lakes_to_plot_22(i)).XY_km_local(:,2), ...
        3000.*ones(length(environs_lakes_2022.boundaries(HF_lakes_to_plot_22(i)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',goldenrod);
end
% MU events
for i=1:1:length(MU_lakes_to_plot_22)
    fill3(environs_lakes_2022.boundaries(MU_lakes_to_plot_22(i)).XY_km_local(:,1), ...
        environs_lakes_2022.boundaries(MU_lakes_to_plot_22(i)).XY_km_local(:,2), ...
         3000.*ones(length(environs_lakes_2022.boundaries(MU_lakes_to_plot_22(i)).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',[0.4 0.8 0.8]);
end
% No-exit freezers 
for i=1:1:length(NE_lakes_to_plot_22)
    fill3(environs_lakes_2022.boundaries(NE_lakes_to_plot_22(i)).XY_km_local(:,1), ...
        environs_lakes_2022.boundaries(NE_lakes_to_plot_22(i)).XY_km_local(:,2), ...
        3000.*ones(length(environs_lakes_2022.boundaries(NE_lakes_to_plot_22(i)).XY_km_local(:,2)),1),...
        'o','MarkerSize',1,'EdgeColor','none','FaceColor',plum(1:3));
end

% HF dates of drainages range (plot as text)
environs_lakes_2022.laketypeing_dates(34,3:4) = [194, 195];
for i=1:1:length(HF_lakes_to_plot_22_INPOLY)
    j=HF_lakes_to_plot_22_INPOLY(i);
   % if ismember(j,HF_lake_gnss_id_2022) % do we have a GNSS drainage date?
   %  [k]=find(HF_lake_gnss_id_2022==j);
   %  text(environs_lakes_2022.X_km(j)-1.0, environs_lakes_2022.Y_km(j),5e3,...
   %      sprintf('%d',HF_lake_gnss_date_2022(k)), ...
   %      'FontSize',fontsize-3,'FontName','Helvetica','FontWeight','bold')
   % else 
       if environs_lakes_2022.laketypeing_dates(j,3)+1 == environs_lakes_2022.laketypeing_dates(j,4)
      text(environs_lakes_2022.X_km(j)-1.0, environs_lakes_2022.Y_km(j),5e3,...
        sprintf('%d',environs_lakes_2022.laketypeing_dates(j,3)+1), ...
        'FontSize',fontsize-3,'FontName','Helvetica','FontWeight','bold')
       else
      text(environs_lakes_2022.X_km(j)-1.5, environs_lakes_2022.Y_km(j),5e3,...
        sprintf('%d\x2013%d',environs_lakes_2022.laketypeing_dates(j,3)+1,environs_lakes_2022.laketypeing_dates(j,4)), ...
        'FontSize',fontsize-3,'FontName','Helvetica','FontWeight','bold')
       %end
   end
end

% Moulin dates of drainages (plot as text)
for i=1:1:length(MU_lakes_to_plot_22_INPOLY)
    j=MU_lakes_to_plot_22_INPOLY(i);
   % if ismember(j,HF_lake_gnss_id_2022) % do we have a GNSS drainage date?
   %  [k]=find(HF_lake_gnss_id_2022==j);
   %  text(environs_lakes_2022.X_km(j)-1.0, environs_lakes_2022.Y_km(j),5e3,...
   %      sprintf('%d',HF_lake_gnss_date_2022(k)), ...
   %      'FontSize',fontsize-3,'FontName','Helvetica','FontWeight','bold')
   % use the S2 catalogue date
   % else 
       if environs_lakes_2022.laketypeing_dates(j,3)+1 == environs_lakes_2022.laketypeing_dates(j,4)
      text(environs_lakes_2022.X_km(j)-1.0, environs_lakes_2022.Y_km(j),5e3,...
        sprintf('%d',environs_lakes_2022.laketypeing_dates(j,3)+1), ...
        'FontSize',fontsize-3,'FontName','Helvetica','FontWeight','bold')
       else 
    text(environs_lakes_2022.X_km(j)-1.5, environs_lakes_2022.Y_km(j),5e3,...
        sprintf('%d\x2013%d',environs_lakes_2022.laketypeing_dates(j,3)+1,environs_lakes_2022.laketypeing_dates(j,4)), ...
        'FontSize',fontsize-3,'FontName','Helvetica','FontWeight','bold')
       end
  % end
end
% 
% % HF ID (plot as text)
% environs_lakes_2022.laketypeing_dates(34,3:4) = [194, 195];
% for i=1:1:length(OS_lakes_to_plot_22)
%     j=OS_lakes_to_plot_22(i);
%    % if ismember(j,HF_lake_gnss_id_2022) % do we have a GNSS drainage date?
%    %  [k]=find(HF_lake_gnss_id_2022==j);
%    %  text(environs_lakes_2022.X_km(j)-1.0, environs_lakes_2022.Y_km(j),5e3,...
%    %      sprintf('%d',HF_lake_gnss_date_2022(k)), ...
%    %      'FontSize',fontsize-3,'FontName','Helvetica','FontWeight','bold')
%    % else 
%        if environs_lakes_2022.laketypeing_dates(j,3)+1 == environs_lakes_2022.laketypeing_dates(j,4)
%       text(environs_lakes_2022.X_km(j), environs_lakes_2022.Y_km(j),5e3,...
%         sprintf('%d',j),'FontSize',fontsize-1,'FontName','Helvetica','FontWeight','bold','Color','r')
%        else
%       text(environs_lakes_2022.X_km(j), environs_lakes_2022.Y_km(j),5e3,...
%         sprintf('%d',j),'FontSize',fontsize-1,'FontName','Helvetica','FontWeight','bold','Color','r')
%        %end
%    end
% end

% site names
% for i=1; text(xy_sta_23_short(i,1)-3.75,xy_sta_23_short(i,2)+0.55,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=2; text(xy_sta_23_short(i,1)-4.25,xy_sta_23_short(i,2)+0.55,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=3; text(xy_sta_23_short(i,1)+1.25,xy_sta_23_short(i,2)+0.55,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=5; text(xy_sta_23_short(i,1)+1.25,xy_sta_23_short(i,2)-0.55,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=6; text(xy_sta_23_short(i,1)+1.25,xy_sta_23_short(i,2)+0.10,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=9; text(xy_sta_23_short(i,1)+1.25,xy_sta_23_short(i,2)+0.55,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=4; text(xy_sta_23_short(i,1)-3.75,xy_sta_23_short(i,2)-0.55,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=7; text(xy_sta_23_short(i,1)-3.75,xy_sta_23_short(i,2)-0.1,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=8; text(xy_sta_23_short(i,1)-3.75,xy_sta_23_short(i,2)+0.75,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=10:11; text(xy_sta_23_short(i,1)-3.75,xy_sta_23_short(i,2)-0.1,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=12:13; text(xy_sta_23_short(i,1)+1.25,xy_sta_23_short(i,2)-0.55,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=14; text(xy_sta_23_short(i,1)-3.75,xy_sta_23_short(i,2)+0.55,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=15:16; text(xy_sta_23_short(i,1)+1.25,xy_sta_23_short(i,2)+0.55,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=17; text(xy_sta_23_short(i,1)+1.25,xy_sta_23_short(i,2)-0.55,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=19:20; text(xy_sta_23_short(i,1)-3.75,xy_sta_23_short(i,2)+0.,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=21:22; text(xy_sta_23_short(i,1)+1.25,xy_sta_23_short(i,2)+0.55,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end

% lake-mechanism example lakes id numbers on the figure
for i=78; text(environs_lakes_2022.X_km(i)-3.0, environs_lakes_2022.Y_km(i)-1.75,5e3,...
        'ID78','FontSize',fontsize-2,'FontName','Helvetica','FontWeight','bold','FontAngle','italic','Color','r'); end
for i=104; text(environs_lakes_2022.X_km(i)+1.35, environs_lakes_2022.Y_km(i)-1.35,5e3,...
        'ID104','FontSize',fontsize-2,'FontName','Helvetica','FontWeight','bold','FontAngle','italic','Color','r'); end
for i=178; text(environs_lakes_2022.X_km(i)+1.25, environs_lakes_2022.Y_km(i)+0.25,5e3,...
        'ID178','FontSize',fontsize-2,'FontName','Helvetica','FontWeight','bold','FontAngle','italic','Color','r'); end
for i=157; text(environs_lakes_2022.X_km(i)-1, environs_lakes_2022.Y_km(i)+1.5,5e3,...
        'ID157','FontSize',fontsize-2,'FontName','Helvetica','FontWeight','bold','FontAngle','italic','Color','r'); end

lgd = legend('Bed Elevation','Ice Elevation','Ice-flow Direction','GNSS station','Moulin',...
    'Overspill drainage','Hydro-fracture drainage','Moulin drainage','No-exit, frozen','Overspill connections');
set(lgd,'Position',[0.805 0.86 .01 0.01],'FontSize',fontsize+1,'FontName','Avenir',...
    'NumColumns',2,'Box','off','TextColor',lake_blue(1:3))

axis equal; grid on;
text(-33, 14.5, 5e3, '2022','FontSize',fontsize+4,'FontName','Helvetica','FontWeight','bold')
set(gca,'FontName','Avenir','FontSize',fontsize,'xtick',-80:10:110,...
    'ytick',-50:10:30,'LineWidth',0.6)
xlim([plot_x(1) plot_x(2)]); ylim([plot_y(1) plot_y(2)]);

axes(axe2) % 2023 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel('East-West [ km ]','FontSize',fontsize);
ylabel('North-South [ km ]','FontSize',fontsize);
set(gca,'xaxislocation','bottom','FontSize',fontsize)

hold on;
% ice-sheet bed contour map, colorbar
surf(BMv5_for_nevis_catchment.X_km/1e3,BMv5_for_nevis_catchment.Y_km./1e3,...
   BMv5_for_nevis_catchment.B_km,'EdgeColor','None','Facealpha',0.1); % plot surface below station points
view([0 90])
clim([-500 500]); colormap(cmap)
text(-40, 14.5, 'c','FontWeight','bold','FontSize',fontsize+3);

% t2=colorbar('east');
% set(t2,'YTick',[-400,-200,0,200,400],'TickDirection','out'); 
% hold all
% set(get(t2,'xlabel'),'String',{'Bed Elevation  [ m a.s.l. ]'},'FontSize',fontsize,'Fontname','Avenir');
% set(t2, 'Position', [0.945 0.885 .006 0.14],'YTick',[-800,-600,-400,-200,0],...
%    'YTickLabel',[-400,-200,0,200,400]);
hold all

% surface ice sheet contours
[C, h] = contour(BMv5_for_nevis_catchment.X_km./1e3, BMv5_for_nevis_catchment.Y_km./1e3, ...
    BMv5_for_nevis_catchment.S_km, surf_contours,'Color',ice_sheet_contours,'LineWidth',0.8);
clabel(C,h,'FontSize',7,'Color',ice_sheet_contours,'FontName','Avenir','LabelSpacing',400)

% ice-flow direction
span = 33; tick_scale = 1.0;
quiver3(BMv5_for_nevis_catchment.X_km(1:span:end,1:span:end)./1e3,...
    BMv5_for_nevis_catchment.Y_km(1:span:end,1:span:end)./1e3,...
    1e3.*ones(size(BMv5_for_nevis_catchment.X_km(1:span:end,1:span:end))),...
    tick_scale.*cosd(flow_angle_interp_2023(1:span:end,1:span:end)),...
    tick_scale.*sind(flow_angle_interp_2023(1:span:end,1:span:end)),...
    0.*sind(flow_angle_interp_2023(1:span:end,1:span:end)),...
    'LineWidth',0.7,'ShowArrowHead','off','Color',[0.5 0.5 0.5],'AutoScale','off');
% [0.6 0.2 0.0] % this is a nice, brown colour, but it's too much

% longterm ELA 
plot(ELA(1074).X_plot, ELA(1074).Y_plot,...
    '-','LineWidth',1.5,'Color',[0.2 0.6 0.4]); % nice green
text(95.5, -33, 10e3, 'ELA','FontName','AvenirBold','FontAngle','Italic',...
    'FontSize',fontsize+1,'Rotation',-90,'Color',[0.2 0.6 0.4]);

% ROI that's relevant to the GNSS array = [-30 100] in X; [-40 10] in Y
% if the lake location is not within relevant catchment, throw it out !! 
plot(ROI_relevant_xv,ROI_relevant_yv,'-','LineWidth',1.5,'Color',[0.6 0.2 0.0]);

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

% count up overspills-->other lakes; overspills-->independent moulins
into_lake_2023 = sum(strcmp(result_23{:,3}, 'L'))
into_independent_moulin_2023 = sum(strcmp(result_23{:,3}, 'M'))

% Natalie Turner's moulins, JUST IN ROI
plot3(moulins_x_2023, moulins_y_2023,300.*ones(length(moulins_x_2023),1),'ko',...
     'MarkerSize',TriangleSize-7,'MarkerFaceColor','w','markerEdgeColor','k',...
     'LineWidth',0.5);
%for i=1:length(moulins_x_2023); text(moulins_x_2023(i), moulins_y_2023(i), 5e3, sprintf('%d',i)); end

% lakes: FASTER S2 2023
[HF_lakes_to_plot] = find(environs_lakes_2023B.laketypeing_dates(:,7)==1); % just HF lakes  
[MU_lakes_to_plot] = find(environs_lakes_2023B.laketypeing_dates(:,7)==2); % just moulin-draining lakes 
[NE_lakes_to_plot] = find(environs_lakes_2023B.laketypeing_dates(:,7)==5); % just no-exit freezers
[OS_lakes_to_plot] = find(environs_lakes_2023B.laketypeing_dates(:,7)==3); % just overspillers
[lakes_to_plot] = vertcat(HF_lakes_to_plot, MU_lakes_to_plot, NE_lakes_to_plot, OS_lakes_to_plot); % all good lakes   

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
for i=1:1:length(NE_lakes_to_plot)
    fill3(environs_lakes_2023B.boundaries(NE_lakes_to_plot(i)).XY_km_local(:,1), ...
        environs_lakes_2023B.boundaries(NE_lakes_to_plot(i)).XY_km_local(:,2), ...
        3000.*ones(length(environs_lakes_2023B.boundaries(NE_lakes_to_plot(i)).XY_km_local(:,2)),1),... 
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',plum(1:3));
end

% lake 3D: 238
fill3(environs_lakes_2023B.boundaries(238).XY_km_local(:,1), ...
        environs_lakes_2023B.boundaries(238).XY_km_local(:,2), ...
        3000.*ones(length(environs_lakes_2023B.boundaries(238).XY_km_local(:,2)),1),... 
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',goldenrod);
% text(environs_lakes_2023B.X_km(238)-1.5, environs_lakes_2023B.Y_km(238),...
%         sprintf('201'), ...
%         'FontSize',fontsize-3,'FontName','Helvetica','FontWeight','bold','Clipping','on')

% Lake 3E: 266+275
% fill3(environs_lakes_2023B_extra.boundaries(266).XY_km_local(:,1), ...
%        environs_lakes_2023B_extra.boundaries(266).XY_km_local(:,2), ...
%         3000.*ones(length(environs_lakes_2023B_extra.boundaries(266).XY_km_local(:,2)),1),...
%         'o','MarkerSize',1,'EdgeColor','none','FaceColor',goldenrod);
fill3(environs_lakes_2023B.boundaries(275).XY_km_local(:,1), ...
        environs_lakes_2023B.boundaries(275).XY_km_local(:,2), ...
         3000.*ones(length(environs_lakes_2023B.boundaries(275).XY_km_local(:,2)),1),...
         'o','MarkerSize',1,'EdgeColor','none','FaceColor',goldenrod);
% text(environs_lakes_2023B.X_km(275)-1.5, environs_lakes_2023B.Y_km(275),...
%         sprintf('201'), ...
%         'FontSize',fontsize-3,'FontName','Helvetica','FontWeight','bold','Clipping','on')

% % L1C from 2022 outline (2022 ID = 77)
% fill(environs_lakes_2022.boundaries(77).XY_km_local(:,1), ...
%         environs_lakes_2022.boundaries(77).XY_km_local(:,2), ...
%          'o','MarkerSize',1,'EdgeColor','none','FaceColor',goldenrod);
% text(environs_lakes_2022.X_km(77)-1.0, environs_lakes_2022.Y_km(77),...
%         sprintf('195'), ...
%         'FontSize',fontsize-3,'FontName','Helvetica','FontWeight','bold','Clipping','on')


% lakes: FASTER S2 2022 -- IN_ROI_POLYGON
clear IN_HF IN_MU IN IN_OS IN_NE
for i=1:length(HF_lakes_to_plot)
    IN_HF(i,1) = inpolygon(environs_lakes_2023B.X_km(HF_lakes_to_plot(i)),environs_lakes_2023B.Y_km(HF_lakes_to_plot(i)),...
        ROI_relevant_xv, ROI_relevant_yv);
end
[HF_lakes_to_plot_INPOLY] = HF_lakes_to_plot(IN_HF);
for i=1:length(MU_lakes_to_plot)
    IN_MU(i,1) = inpolygon(environs_lakes_2023B.X_km(MU_lakes_to_plot(i)),environs_lakes_2023B.Y_km(MU_lakes_to_plot(i)),...
        ROI_relevant_xv, ROI_relevant_yv);
end
[MU_lakes_to_plot_INPOLY] = MU_lakes_to_plot(IN_MU);
for i=1:length(NE_lakes_to_plot)
    IN_NE(i,1) = inpolygon(environs_lakes_2023B.X_km(NE_lakes_to_plot(i)),environs_lakes_2023B.Y_km(NE_lakes_to_plot(i)),...
        ROI_relevant_xv, ROI_relevant_yv);
end
[NE_lakes_to_plot_INPOLY] = NE_lakes_to_plot(IN_NE);
for i=1:length(lakes_to_plot)
    IN(i,1) = inpolygon(environs_lakes_2023B.X_km(lakes_to_plot(i)),environs_lakes_2023B.Y_km(lakes_to_plot(i)),...
        ROI_relevant_xv, ROI_relevant_yv);
end
[lakes_to_plot_23_INPOLY] = lakes_to_plot(IN);
for i=1:length(OS_lakes_to_plot)
    IN_OS(i,1) = inpolygon(environs_lakes_2023B.X_km(OS_lakes_to_plot(i)),environs_lakes_2023B.Y_km(OS_lakes_to_plot(i)),...
        ROI_relevant_xv, ROI_relevant_yv);
end
[OS_lakes_to_plot_23_INPOLY] = OS_lakes_to_plot(IN_OS);
% for mapping overspill fates with ROI:
mask_overspill = false(406,1); mask_overspill(OS_lakes_to_plot_23_INPOLY) = true;

% lake stats: count up overspills-->other lakes; overspills-->terminal moulins
ROI_area_23 = polyarea(ROI_relevant_xv,ROI_relevant_yv)
in_ROI_lakes_2023 = length(lakes_to_plot_23_INPOLY)
lake_INPOLY_density_23 = length(lakes_to_plot_23_INPOLY)./ROI_area
HF_in_ROI_lakes_2023 = length(HF_lakes_to_plot_INPOLY)
MU_in_ROI_lakes_2023 = length(MU_lakes_to_plot_INPOLY)
OS_in_ROI_lakes_2023 = length(OS_lakes_to_plot_23_INPOLY)
NE_in_ROI_lakes_2023 = length(NE_lakes_to_plot_INPOLY)
OS_into_lake_2023 = sum(strcmp(result_23{OS_lakes_to_plot_23_INPOLY,3}, 'L'))
OS_into_independent_moulin_2023 = sum(strcmp(result_23{OS_lakes_to_plot_23_INPOLY,3}, 'M'))

% HF dates of drainages range (plot as text)
for i=1:1:length(HF_lakes_to_plot_INPOLY)
    j=HF_lakes_to_plot_INPOLY(i);
   % if ismember(j,HF_lake_gnss_id_2023) % do we have a GNSS drainage date?
   %  [k]=find(HF_lake_gnss_id_2023==j);
   %  text(environs_lakes_2023.X_km_origin(j)-1.0, environs_lakes_2023.Y_km_origin(j),5e3,...
   %      sprintf('%d',HF_lake_gnss_date_2023(k)), ...
   %      'FontSize',fontsize-3,'FontName','Helvetica','FontWeight','bold','Clipping','on')
   % else 
       if environs_lakes_2023B.laketypeing_dates(j,3)+1 == environs_lakes_2023B.laketypeing_dates(j,4)
      text(environs_lakes_2023B.X_km(j)-1.0, environs_lakes_2023B.Y_km(j),5e3,...
        sprintf('%d',environs_lakes_2023B.laketypeing_dates(j,3)+1), ...
        'FontSize',fontsize-3,'FontName','Helvetica','FontWeight','bold')
    else
      text(environs_lakes_2023B.X_km(j)-1.5, environs_lakes_2023B.Y_km(j),5e3,...
        sprintf('%d\x2013%d',environs_lakes_2023B.laketypeing_dates(j,3)+1,environs_lakes_2023B.laketypeing_dates(j,4)), ...
        'FontSize',fontsize-3,'FontName','Helvetica','FontWeight','bold')
       end
   % end
end

% Moulin dates of drainages (plot as text)
for i=1:1:length(MU_lakes_to_plot_INPOLY)
    j=MU_lakes_to_plot_INPOLY(i);
   % if ismember(j,HF_lake_gnss_id_2023) % do we have a GNSS drainage date?
   %  [k]=find(HF_lake_gnss_id_2023==j);
   %  text(environs_lakes_2023.X_km(j)-1.0, environs_lakes_2023.Y_km(j),5e3,...
   %      sprintf('%d',HF_lake_gnss_date_2023(k)), ...
   %      'FontSize',fontsize-3,'FontName','Helvetica','FontWeight','bold')
   % % use the S2 catalogue date
   % else 
       if environs_lakes_2023B.laketypeing_dates(j,3)+1 == environs_lakes_2023B.laketypeing_dates(j,4)
      text(environs_lakes_2023B.X_km(j)-1.0, environs_lakes_2023B.Y_km(j),5e3,...
        sprintf('%d',environs_lakes_2023B.laketypeing_dates(j,3)+1), ...
        'FontSize',fontsize-3,'FontName','Helvetica','FontWeight','bold')
        else 
    text(environs_lakes_2023B.X_km(j)-1.5, environs_lakes_2023B.Y_km(j),5e3,...
        sprintf('%d\x2013%d',environs_lakes_2023B.laketypeing_dates(j,3)+1,environs_lakes_2023B.laketypeing_dates(j,4)), ...
        'FontSize',fontsize-3,'FontName','Helvetica','FontWeight','bold')
        end
   % end
end

% main array GNSS
plot3(xy_sta_23_short(:,1),xy_sta_23_short(:,2),700.*ones(num_sta,1),'k^',...
    'MarkerSize',TriangleSize-3,'MarkerFaceColor','w','markerEdgeColor','k',...
    'LineWidth',0.8); 

% site names
% for i=1; text(xy_sta_23_short(i,1)-3.75,xy_sta_23_short(i,2)+0.55,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=2; text(xy_sta_23_short(i,1)-4.25,xy_sta_23_short(i,2)+0.55,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=3; text(xy_sta_23_short(i,1)+1.25,xy_sta_23_short(i,2)+0.55,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=5; text(xy_sta_23_short(i,1)+1.25,xy_sta_23_short(i,2)-0.55,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=6; text(xy_sta_23_short(i,1)+1.25,xy_sta_23_short(i,2)+0.10,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=9; text(xy_sta_23_short(i,1)+1.25,xy_sta_23_short(i,2)+0.55,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=4; text(xy_sta_23_short(i,1)-3.75,xy_sta_23_short(i,2)-0.55,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=7; text(xy_sta_23_short(i,1)-3.75,xy_sta_23_short(i,2)-0.1,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=8; text(xy_sta_23_short(i,1)-3.75,xy_sta_23_short(i,2)+0.75,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=10:11; text(xy_sta_23_short(i,1)-3.75,xy_sta_23_short(i,2)-0.1,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=12:13; text(xy_sta_23_short(i,1)+1.25,xy_sta_23_short(i,2)-0.55,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=14; text(xy_sta_23_short(i,1)-3.75,xy_sta_23_short(i,2)+0.55,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=15:16; text(xy_sta_23_short(i,1)+1.25,xy_sta_23_short(i,2)+0.55,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=17:18; text(xy_sta_23_short(i,1)+1.25,xy_sta_23_short(i,2)-0.55,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=19:20; text(xy_sta_23_short(i,1)-3.75,xy_sta_23_short(i,2)+0.,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end
% for i=21:22; text(xy_sta_23_short(i,1)+1.25,xy_sta_23_short(i,2)+0.55,station_names_short(i),'Fontname','AvenirBold','FontSize',fontsize-2,'FontAngle','italic'); end

axis equal; grid on; box on;
text(-33, 14.5, 5e3, '2023','FontSize',fontsize+4,'FontName','Helvetica','FontWeight','bold')
xlim([plot_x(1) plot_x(2)]); ylim([plot_y(1) plot_y(2)]);
set(gca,'FontName','Avenir','FontSize',fontsize,'xtick',-80:10:110,...
    'ytick',-50:10:30,'LineWidth',0.6)

% %% mini Green Greenland map
% axe_map = axes('Position',[0.845 0.840 0.17 0.17],'Box','off',...
%     'XTickLabel',[],'YTickLabel',[]);
% axes(axe_map)
% ax = worldmap('greenland');
% greenland = shaperead('landareas', 'UseGeoCoords', true,...
%   'Selector',{@(name) strcmp(name,'Greenland'), 'Name'});
% patchm(greenland.Lat,greenland.Lon,'FaceColor', [0.4 0.6 0.4],'EdgeColor',lake_blue(1:3),...
%     'LineWidth',1.1); hold on;
% %plotm([69],[-49],'p','MarkerSize',20,'MarkerFaceColor',lake_blue(1:3),'MarkerEdgeColor',[1 1 1]); 
% fillm([68.2 69.2 69.2 68.2],[-47.5 -47.5 -50.5 -50.5],'w','LineWidth',1.1,'FaceColor',goldenrod,'EdgeColor',goldenrod)
% axis off; framem off; gridm off; mlabel off; plabel off; 

%% print figure
print(gcf,'-dpng','-r300',sprintf('../../paperfigs/paperfig1_map_fullROI_260130.png')); 

