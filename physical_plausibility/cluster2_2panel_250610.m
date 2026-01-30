% physics-informed cluster chronology exploration:
% map-view-a: lakes that HF in cluster (lake margins goldenrod)
%            lakes that are voluminous enough to HF
%            subglacial q for first possible day of cluster (colormap)
%            lakes that have already drained (grey)
% map-view-b: lakes that HF in cluster (lake margins goldenrod)
%            lakes that are voluminous enough to HF 
%            sigma1 from elastic blister opening + background flow (colormap)
%            sigma1 theshold contour (200 kPa; black)
%            rings moving out from HF drainages at flood-propagation speed
%            lakes that have already drained (grey)

% 24-09-18 LAS: physics-informed cluster chronology first try for C1 (2022)
% 24-09-18 LAS: add in elastic blister stress 
% 24-09-24 LAS: add in flow directions
% 24-10-14 LAS: this is C2 now, let's be chronological in cluster naming
% 24-10-20 LAS: try pixel count of max surface area --> lake volume
% 25-02-19 LAS: correct indexing for \sigma_{1,bg}
% 25-03-02 LAS: switch to one, three-panel figure
% 25-04-23 LAS: update to polarstereographic lake locations
% 25-06-10 LAS: remove scatter, sigma_{1,elastic} only

clear all; close all;

%% load in 2022 preliminary datasets needed for cluster chrono
preliminaries_2022
num_sta_22 = num_sta_22-2;
% SET DAY OF NEVIS TO BE PLOTTED (as first day of HF window)
DOY = 210; % DOY

%% ELASTIC STRESS CHANGE DUE TO BLISTER
addpath('make_idealized_blisters/','forward_okada85/');
% lake volumes for C2 cluster:
cluster_name = 'C2';
cluster_mid = [168, 173, 183, 185, 201, 223];

% lake volume: option 1 -- linear model of lake surface-area growth:
% volume_mid = environs_lakes_2022.lake_volume_interp(cluster_mid,DOY-1); % m^3 [ linear model]

% lake volume: option 2 -- use maximum surface (pixel) area of the lake:
% assume conical lake; assign an average diameter to the lake
lake_diameter_interp = 2.*sqrt((environs_lakes_2022.max_surface_area.*1e6)./3.14159); % average diameter of this lake [m]
% crude volume estimation (pi * r^2 * (h/3)), assuming lake-bottom slope 
% that's set by the ratio:  [ diameter:cone_height :: 100:1 ]
lake_volume_interp = 3.14159.*((lake_diameter_interp./2).^2).* ...
        ((lake_diameter_interp./100)./3);  % [m^3]
volume_mid = lake_volume_interp(cluster_mid); % m^3 [ max pixel count of lake ]

ice_thickness_mid = environs_lakes_2022.H(cluster_mid); % m
% all thicknesses are within +/-75 m of 1250m, so use 1.25-km for blister depth 
% blister radius, height, slip for these lakes:
[patches_idealized] = idealized_blister_function(volume_mid,cluster_name); % input volumes in m^3 !!!
% forward Okada [1985] for these lakes, taking H=1.25km:
[sigma_ij_1250m_SLder] = Nsurface_displ_strain_stress_1250m_function(patches_idealized,cluster_name);

%% create map of ROI with these elastic stresses in the right place
% origins for each lake:
origin_cluster = [environs_lakes_2022.X_km(cluster_mid); environs_lakes_2022.Y_km(cluster_mid)];
% blank X-Y grid grid points = nevis xx yy
xx_vec = reshape(xx,601*264,1); 
yy_vec = reshape(yy,601*264,1);
for i=1:length(cluster_mid)
    % flow angle
    theta(i,1) = flow_angle((ID_stresses_sig1_winter(cluster_mid(i),1))); % [ deg ]
    % rotation matrix
    Rotation=[cosd(theta(i,1)) -sind(theta(i,1)); sind(theta(i,1)) cosd(theta(i,1))]; % create rotation matrix
    % align with ice-flow direction by rotating grid
    rotXY = Rotation*[sigma_ij_1250m_SLder.xx(:).'; sigma_ij_1250m_SLder.yy(:).']; % vectors X rotation matrix 
    XX = reshape(rotXY(1, :), size(sigma_ij_1250m_SLder.xx));
    YY = reshape(rotXY(2, :), size(sigma_ij_1250m_SLder.yy));
    % translate tau_xx output (x_star) --> origin of lake location in nevis xx, yy
    x_sigma = reshape(XX,121*121,1) + origin_cluster(1,i); % [ km ]
    y_sigma = reshape(YY,121*121,1) + origin_cluster(2,i); % [ km ]
    % interpolate stress output to nevis grid
    nevis_sigma_xx(i,:) = griddata(x_sigma,y_sigma,sigma_ij_1250m_SLder.sigma_xx(i,:),xx_vec,yy_vec); % xx
    nevis_sigma_yy(i,:) = griddata(x_sigma,y_sigma,sigma_ij_1250m_SLder.sigma_yy(i,:),xx_vec,yy_vec); % yy
    nevis_sigma_xy(i,:) = griddata(x_sigma,y_sigma,sigma_ij_1250m_SLder.sigma_xy(i,:),xx_vec,yy_vec); % xy
end
% add all lakes' elastic stress together:
nevis_sigma_xx_cluster = nansum(nevis_sigma_xx); % Pa
nevis_sigma_yy_cluster = nansum(nevis_sigma_yy); % Pa
nevis_sigma_xy_cluster = nansum(nevis_sigma_xy); % Pa
% calculate sig1_elastic:
princ_sigma1_elastic = (0.5.*(nevis_sigma_xx_cluster+nevis_sigma_yy_cluster)) + ...
 (sqrt(((0.5.*(nevis_sigma_xx_cluster-nevis_sigma_yy_cluster)).^2) + (nevis_sigma_xy_cluster.^2))); 

v200=200;
% % check figs :D
figure;
contourf(xx,yy,reshape(princ_sigma1_elastic./1e3,601,264)); hold on
contour(xx,yy,reshape(princ_sigma1_elastic./1e3,601,264),[v200 v200],'k','LineWidth',1.2);
colorbar; shg

for i=4
princ_sigma1_elastic_1 = (0.5.*(nevis_sigma_xx(i,:)+nevis_sigma_yy(i,:))) + ...
 (sqrt(((0.5.*(nevis_sigma_xx(i,:)-nevis_sigma_yy(i,:))).^2) + (nevis_sigma_xy(i,:).^2))); 
end

figure;
contourf(xx,yy,reshape(princ_sigma1_elastic_1./1e3,601,264)); hold on
contour(xx,yy,reshape(princ_sigma1_elastic_1./1e3,601,264),[v200 v200],'k','LineWidth',1.2);
colorbar; axis equal; shg

figure;
contourf(xx,yy,reshape(nevis_sigma_yy(i,:)./1e3,601,264)); hold on
contour(xx,yy,reshape(nevis_sigma_yy(i,:)./1e3,601,264),[v200 v200],'k','LineWidth',1.2);
colorbar; axis equal; shg

%% BLISTER PROPAGATION DIRECTION WITH TIME
% blister_speed = 0.3; % m/s [hoffman [2018] rough average for WGrIS]
blister_speed = 0.4; % m/s [ a good average for our three observed floods ]
% What direction is down subglacial flood path?
clear theta_range distance_propagated time_range x_circle y_circle
theta_range = 140:1:210; % degrees with East = 0 degrees. --> 230 to“310 in ice-flow direction
% time of calculations [12 hr, 24 hr, 36 hr, 48 hr]
time_range = (12:12:48).*(60*60); % [s]
% distance travelled = radius
distance_propagated = (blister_speed.*time_range)./1e3; % km 
% circle origins = HF lake origins :)
for i=1:length(cluster_mid); for j=1:length(time_range)
    % circle xy-points
    x_circle(i,j,:) = origin_cluster(1,i) + (distance_propagated(j).*cosd(theta_range));
    y_circle(i,j,:) = origin_cluster(2,i) + (distance_propagated(j).*sind(theta_range));
end; end
% test figure
figure; clf; 
for i=1:length(cluster_mid); for j=1:length(time_range)
plot(squeeze(x_circle(i,j,:)),squeeze(y_circle(i,j,:)),'-'); hold on;
end; end

%% CLUSTER CHRONO MAP FIG
close all
% figure preliminaries
fontsize = 10-1;  mm=fontsize;
TriangleSize = 5;
DiamondSize = 1.5;
HF_size = 70;
lake_size = 30;
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
% more lake colors
filling_color = [0.2 0.6 1.0];
draining_color = [1.0 0.8 0.0];
still_draining_color = [0.2 0.8 0.2];
frozen_color = [0.8 0.8 0.8];

% lake drainage dates constrained by GNSS
HF_lake_gnss_id_2022 = [78, 81, 143, 156, 226, 168]; % 1A, 1B, 2A, 2B, 3B, MHIH-lake
HF_lake_gnss_date_2022 = [195, 195, 214, 214, 209, 209]; % 1A, 1B, 2A, 2B, 3B, MHIH-lake
HF_lake_gnss_id_2023 = [81, 86, 132, 142, 201, 208, 150]; % 1A, 1B, 2A, 2B, 3A, 3B, MHIH-lake
HF_lake_gnss_date_2023 = [191, 191, 202, 202, 200, 200, 196]; % 1A, 1B, 2A, 2B, 3A, 3B, MHIH-lake

figure(1); clf; 
set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',1.5.*[0.5 0 15.3 15]);
axeA = axes('Position',[0.048 0.53 0.91 0.46],'Box','on','NextPlot','add','XTickLabel',[]);
axeB = axes('Position',[0.048 0.05 0.91 0.46],'Box','on','NextPlot','add');
% full ROI best for subglacial discharge
plot_x = [-45 90];
plot_y = [-50 17];

axes(axeA) % 2022 subglacial discharge %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ylabel('North-South [ km ]','FontSize',fontsize);
set(gca,'xaxislocation','bottom','FontSize',fontsize)
text(-51, 16, 'a','FontWeight','bold','FontSize',fontsize+3);

hold all

% nevis discharge: plot frames
    disp(['Frame ',num2str(DOY),' / ',num2str(length(DOY)),' ...']);
    % load timestep
    load([fn,'/',int2four(DOY)]);
    % extract new variables
    if isfield(vv,'nbdy'), gg = nevis_label(gg,vv.nbdy); aa.phi = aa.phi_a(gg.nbdy); end % remask
    aa = nevis_inputs(vv.t,aa,pp,gg,oo);
    oo.evaluate_variables = 1; oo.evaluate_residual = 0; oo.evaluate_jacobian = 0; [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
    vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
    nevis_unpack(aa,vv2);
    clear tt vv vv2

    % discharge plotting
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny;
    zz = ps.qs*reshape(qs+qe+qQ,gg.nI,gg.nJ);   
    cax = [10^(-4) 10^(0)]; 
    zz(zz<10^(-4)) = 0; % don't plot low-flow areas 
    % logarithmic pcolor
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis
    load cmapq2; cmap2=cmap; colormap(axeA, cmap2); %cmap2(64,1:3) = NaN;
    hand = pcolor(xx,yy,log10(zz)); % shading interp;
    set(hand,'linestyle','none','FaceAlpha',1.0); 
    caxis(log10(cax)); 
    axis image;
    % color bar
    cbar_size = [0.77 0.546 0.165 0.005];
    cx5 = colorbar('peer',gca,'horizontal','position',cbar_size,'XAxislocation','top');
    tmp = get(cx5,'XTick');  labels = 10.^tmp;  
    xlabel(cx5, 'Subglacial Discharge {\itq} [ m^2 s^{-1} ]')
    set(cx5,'XTickLabel',labels,'tickdir','out');

% surface ice sheet contours
surf_contours = 0:100:1900; % m a.s.l.
[C, h] = contour(BMv5_for_nevis_catchment_noSK.X_m./1e3, BMv5_for_nevis_catchment_noSK.Y_m./1e3, ...
    BMv5_for_nevis_catchment_noSK.S_m, surf_contours,'Color',ice_sheet_contours,'LineWidth',0.8);
clabel(C,h,'FontSize',7,'Color',ice_sheet_contours,'FontName','Avenir','LabelSpacing',400)

% % ROI that's relevant to the GNSS array = [-30 100] in X; [-40 10] in Y
% % if the lake location is not within relevant catchment, throw it out !! 
plot(ROI_relevant_xv,ROI_relevant_yv,'-','LineWidth',2,'Color',[0.6 0.2 0.0]);

% main array GNSS
plot3(xy_sta_22_short(:,1),xy_sta_22_short(:,2),300.*ones(num_sta_22,1),'k^',...
    'MarkerSize',4,'MarkerFaceColor','k','markerEdgeColor','k');

% cluster cOneXiOnS
grey_connections = [0 0 0];
% mid C3
cluster_mid = [168, 173, 183, 185, 201, 223];
cc3 = plot([environs_lakes_2022.X_km(cluster_mid(1)) environs_lakes_2022.X_km(cluster_mid(2))],...
    [environs_lakes_2022.Y_km(cluster_mid(1)) environs_lakes_2022.Y_km(cluster_mid(2))],...
    '-','LineWidth',2,'Color',grey_connections);
cc31 =plot([environs_lakes_2022.X_km(cluster_mid(2)) environs_lakes_2022.X_km(cluster_mid(3))],...
    [environs_lakes_2022.Y_km(cluster_mid(2)) environs_lakes_2022.Y_km(cluster_mid(3))],...
    '-','LineWidth',2,'Color',grey_connections);
cc32 =plot([environs_lakes_2022.X_km(cluster_mid(4)) environs_lakes_2022.X_km(cluster_mid(5))],...
    [environs_lakes_2022.Y_km(cluster_mid(4)) environs_lakes_2022.Y_km(cluster_mid(5))],...
    '-','LineWidth',2,'Color',grey_connections);
cc33 =plot([environs_lakes_2022.X_km(cluster_mid(6)) environs_lakes_2022.X_km(cluster_mid(5))],...
    [environs_lakes_2022.Y_km(cluster_mid(6)) environs_lakes_2022.Y_km(cluster_mid(5))],...
    '-','LineWidth',2,'Color',grey_connections);
% lake identification number labels
text(environs_lakes_2022.X_km(cluster_mid(1))+1.25, environs_lakes_2022.Y_km(cluster_mid(1))-0.5,...
    "L168","Color",'k',"FontSize",fontsize,"FontName","Helvetical",'FontWeight','Bold','FontAngle','italic')
text(environs_lakes_2022.X_km(cluster_mid(2))+1, environs_lakes_2022.Y_km(cluster_mid(2)),...
    "L173","Color",'k',"FontSize",fontsize,"FontName","Helvetical",'FontWeight','Bold','FontAngle','italic')
text(environs_lakes_2022.X_km(cluster_mid(3))+1, environs_lakes_2022.Y_km(cluster_mid(3)),...
    "L183","Color",'k',"FontSize",fontsize,"FontName","Helvetical",'FontWeight','Bold','FontAngle','italic')
text(environs_lakes_2022.X_km(cluster_mid(4))+0.25, environs_lakes_2022.Y_km(cluster_mid(4))+2,...
    "L185","Color",'k',"FontSize",fontsize,"FontName","Helvetical",'FontWeight','Bold','FontAngle','italic')
text(environs_lakes_2022.X_km(cluster_mid(5))+1, environs_lakes_2022.Y_km(cluster_mid(5))+1.25,...
    "L201","Color",'k',"FontSize",fontsize,"FontName","Helvetical",'FontWeight','Bold','FontAngle','italic')
text(environs_lakes_2022.X_km(cluster_mid(6))+1, environs_lakes_2022.Y_km(cluster_mid(6))-0.25,...
    "L223","Color",'k',"FontSize",fontsize,"FontName","Helvetical",'FontWeight','Bold','FontAngle','italic')

% CLASSIFIER PLOTTING
% lake filling
Filling = daily_lake_classifier(:,DOY)==1;
% HF lake filling
scatter(X_km(1,type_num==1 & Filling==1),Y_km(1,type_num==1 & Filling==1),HF_size,...
    'p','filled','MarkerFaceColor',filling_color,'MarkerEdgeColor','k','LineWidth',0.8); % HF
% moulin lake filling
scatter(X_km(1,type_num==2 & Filling==1),Y_km(1,type_num==2 & Filling==1),lake_size,...
    'v','filled','MarkerFaceColor',filling_color,'MarkerEdgeColor','k'); % moulin
% overspill lake filling
scatter(X_km(1,type_num==3 & Filling==1),Y_km(1,type_num==3 & Filling==1),lake_size,...
    'o','filled','MarkerFaceColor',filling_color,'MarkerEdgeColor','k'); % overspill
% crevasses filling
scatter(X_km(1,type_num==4 & Filling==1),Y_km(1,type_num==4 & Filling==1),lake_size,...
    'd','filled','MarkerFaceColor',filling_color,'MarkerEdgeColor','k'); % crevasses
% frozen lake filling
scatter(X_km(1,type_num==5 & Filling==1),Y_km(1,type_num==5 & Filling==1),lake_size,...
    's','filled','MarkerFaceColor',filling_color,'MarkerEdgeColor','k'); % frozen

% lake draining -- initial
Draining = daily_lake_classifier(:,DOY)==2;
% HF lake draining
scatter3(X_km(1,type_num==1 & Draining==1),Y_km(1,type_num==1 & Draining==1),...
    1e10.*abs(X_km(1,type_num==1 & Draining==1)),HF_size+40,...
    'p','filled','MarkerFaceColor',draining_color,'MarkerEdgeColor','k'); % HF
% moulin lake draining
scatter(X_km(1,type_num==2 & Draining==1),Y_km(1,type_num==2 & Draining==1),lake_size,...
    'v','filled','MarkerFaceColor',draining_color,'MarkerEdgeColor','k'); % moulin
% overspill lake draining
scatter(X_km(1,type_num==3 & Draining==1),Y_km(1,type_num==3 & Draining==1),lake_size,...
    'o','filled','MarkerFaceColor',draining_color,'MarkerEdgeColor','k'); % overspill
% crevasses draining
scatter(X_km(1,type_num==4 & Draining==1),Y_km(1,type_num==4 & Draining==1),lake_size,...
    'd','filled','MarkerFaceColor',draining_color,'MarkerEdgeColor','k'); % crevasses

% lake draining -- continued
Still_Draining = daily_lake_classifier(:,DOY)==3;
% moulin lake draining
scatter(X_km(1,type_num==2 & Still_Draining==1),Y_km(1,type_num==2 & Still_Draining==1),lake_size,...
    'v','filled','MarkerFaceColor',still_draining_color,'MarkerEdgeColor','k'); % moulin
% overspill lake draining
scatter(X_km(1,type_num==3 & Still_Draining==1),Y_km(1,type_num==3 & Still_Draining==1),lake_size,...
    'o','filled','MarkerFaceColor',still_draining_color,'MarkerEdgeColor','k'); % overspill
% crevasses draining
scatter(X_km(1,type_num==4 & Still_Draining==1),Y_km(1,type_num==4 & Still_Draining==1),lake_size,...
    'd','filled','MarkerFaceColor',still_draining_color,'MarkerEdgeColor','k'); % crevasses

% HF empty basin
Empty = daily_lake_classifier(:,DOY)==3;
% HF empty basin
scatter3(X_km(1,type_num==1 & Empty==1),Y_km(1,type_num==1 & Empty==1),...
    1e10.*abs(X_km(1,type_num==1 & Empty==1)),HF_size,...
    'p','filled','MarkerFaceColor',frozen_color,'MarkerEdgeColor','k','LineWidth',0.8); % HF

% lake frozen
Frozen = daily_lake_classifier(:,DOY)==4;
% moulin lake frozen
scatter(X_km(1,type_num==2 & Frozen==1),Y_km(1,type_num==2 & Frozen==1),lake_size,...
    'v','filled','MarkerFaceColor',frozen_color,'MarkerEdgeColor','k','LineWidth',0.8); % moulin
% overspill lake frozen
scatter(X_km(1,type_num==3 & Frozen==1),Y_km(1,type_num==3 & Frozen==1),lake_size,...
    'o','filled','MarkerFaceColor',frozen_color,'MarkerEdgeColor','k','LineWidth',0.8); % overspill
% crevasses frozen
scatter(X_km(1,type_num==4 & Frozen==1),Y_km(1,type_num==4 & Frozen==1),lake_size,...
    'd','filled','MarkerFaceColor',frozen_color,'MarkerEdgeColor','k','LineWidth',0.8); % crevasses
% frozen lake frozen
scatter(X_km(1,type_num==5 & Frozen==1),Y_km(1,type_num==5 & Frozen==1),lake_size,...
    's','filled','MarkerFaceColor',frozen_color,'MarkerEdgeColor','k','LineWidth',0.8); % frozen

% HOW MANY VIABLE LAKES?
count_num_lakes_viable_to_HF = sum(Filling) + sum(Draining) + sum(Still_Draining)

lgdA = legend('Subglacial Discharge {\itq}','Ice Elevation','Subglacial ROI',...
    'GNSS station','C2 HF events');
set(lgdA,'Position',[0.885 0.930 0.05 0.03],'FontSize',fontsize-3,'FontName','Avenir',...
    'NumColumns',1,'Box','on','TextColor',lake_blue(1:3))

axis equal; grid on;
date_of_q = sprintf('2022/%d',DOY);
text(-44, 15, date_of_q,'FontSize',fontsize+4,'FontName','Helvetica','FontWeight','bold')
set(gca,'FontName','Avenir','FontSize',fontsize,'xtick',-80:10:110,...
    'ytick',-50:10:30,'LineWidth',0.6,'Layer','top')
xlim([plot_x(1) plot_x(2)]); ylim([plot_y(1) plot_y(2)]);

% sigma1 map ELASTIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(axeB)
xlabel('East-West [ km ]','FontSize',fontsize);
ylabel('North-South [ km ]','FontSize',fontsize);
set(gca,'xaxislocation','bottom','FontSize',fontsize)
text(-51, 16, 'b','FontWeight','bold','FontSize',fontsize+3);

hold all
% elastic blister sigma1 + bound viscous
levels_kPa=(-1800:10:1800);
[~,H] = contourf(xx,yy,(reshape(princ_sigma1_elastic./1e3,601,264)),levels_kPa);
set(H, 'LineColor','none');
[~,H2] = contour(xx,yy,(reshape(princ_sigma1_elastic./1e3,601,264)),[200 200],'k','LineWidth',1.2);
set(H2, 'LineColor','k');
% color bar
    caxis([-600 600])
    colormap(axeB,bluewhitered(100));
    cbar_size = [0.770 0.062 0.165 0.005];
    cx1 = colorbar('peer',gca,'horizontal','position',cbar_size,'XAxislocation','top');
    xlabel(cx1, '\sigma_{1}  [ kPa ]')
    set(cx1,'XTick',[-600 -300 0 300 600],'tickdir','out');

% propagating blister speed radar circles:
i=1; plot(squeeze(x_circle(i,2,:)),squeeze(y_circle(i,2,:)),'-','LineWidth',2,'Color',[0.6 0.6 0.6]); 

% CLASSIFIER PLOTTING
% lake filling
Filling = daily_lake_classifier(:,DOY+2)==1;
% HF lake filling
scatter(X_km(1,type_num==1 & Filling==1),Y_km(1,type_num==1 & Filling==1),HF_size,...
    'p','filled','MarkerFaceColor',filling_color,'MarkerEdgeColor','k','LineWidth',0.8); % HF
% moulin lake filling
scatter(X_km(1,type_num==2 & Filling==1),Y_km(1,type_num==2 & Filling==1),lake_size,...
    'v','filled','MarkerFaceColor',filling_color,'MarkerEdgeColor','k'); % moulin
% overspill lake filling
scatter(X_km(1,type_num==3 & Filling==1),Y_km(1,type_num==3 & Filling==1),lake_size,...
    'o','filled','MarkerFaceColor',filling_color,'MarkerEdgeColor','k'); % overspill
% crevasses filling
scatter(X_km(1,type_num==4 & Filling==1),Y_km(1,type_num==4 & Filling==1),lake_size,...
    'd','filled','MarkerFaceColor',filling_color,'MarkerEdgeColor','k'); % crevasses
% frozen lake filling
scatter(X_km(1,type_num==5 & Filling==1),Y_km(1,type_num==5 & Filling==1),lake_size,...
    's','filled','MarkerFaceColor',filling_color,'MarkerEdgeColor','k'); % frozen

% lake draining -- initial
Draining = daily_lake_classifier(:,DOY+2)==2;
% HF lake draining
scatter3(X_km(1,type_num==1 & Draining==1),Y_km(1,type_num==1 & Draining==1),...
    1e10.*abs(X_km(1,type_num==1 & Draining==1)),HF_size+40,...
    'p','filled','MarkerFaceColor',draining_color,'MarkerEdgeColor','k'); % HF
% moulin lake draining
scatter(X_km(1,type_num==2 & Draining==1),Y_km(1,type_num==2 & Draining==1),lake_size,...
    'v','filled','MarkerFaceColor',draining_color,'MarkerEdgeColor','k'); % moulin
% overspill lake draining
scatter(X_km(1,type_num==3 & Draining==1),Y_km(1,type_num==3 & Draining==1),lake_size,...
    'o','filled','MarkerFaceColor',draining_color,'MarkerEdgeColor','k'); % overspill
% crevasses draining
scatter(X_km(1,type_num==4 & Draining==1),Y_km(1,type_num==4 & Draining==1),lake_size,...
    'd','filled','MarkerFaceColor',draining_color,'MarkerEdgeColor','k'); % crevasses

% lake draining -- continued
Still_Draining = daily_lake_classifier(:,DOY+2)==3;
% moulin lake draining
scatter(X_km(1,type_num==2 & Still_Draining==1),Y_km(1,type_num==2 & Still_Draining==1),lake_size,...
    'v','filled','MarkerFaceColor',still_draining_color,'MarkerEdgeColor','k'); % moulin
% overspill lake draining
scatter(X_km(1,type_num==3 & Still_Draining==1),Y_km(1,type_num==3 & Still_Draining==1),lake_size,...
    'o','filled','MarkerFaceColor',still_draining_color,'MarkerEdgeColor','k'); % overspill
% crevasses draining
scatter(X_km(1,type_num==4 & Still_Draining==1),Y_km(1,type_num==4 & Still_Draining==1),lake_size,...
    'd','filled','MarkerFaceColor',still_draining_color,'MarkerEdgeColor','k'); % crevasses

% HF empty basin
Empty = daily_lake_classifier(:,DOY+2)==3;
% HF empty basin
scatter3(X_km(1,type_num==1 & Empty==1),Y_km(1,type_num==1 & Empty==1),...
    1e10.*abs(X_km(1,type_num==1 & Empty==1)),HF_size,...
    'p','filled','MarkerFaceColor',frozen_color,'MarkerEdgeColor','k','LineWidth',0.8); % HF

% lake frozen
Frozen = daily_lake_classifier(:,DOY+2)==4;
% moulin lake frozen
scatter(X_km(1,type_num==2 & Frozen==1),Y_km(1,type_num==2 & Frozen==1),lake_size,...
    'v','filled','MarkerFaceColor',frozen_color,'MarkerEdgeColor','k','LineWidth',0.8); % moulin
% overspill lake frozen
scatter(X_km(1,type_num==3 & Frozen==1),Y_km(1,type_num==3 & Frozen==1),lake_size,...
    'o','filled','MarkerFaceColor',frozen_color,'MarkerEdgeColor','k','LineWidth',0.8); % overspill
% crevasses frozen
scatter(X_km(1,type_num==4 & Frozen==1),Y_km(1,type_num==4 & Frozen==1),lake_size,...
    'd','filled','MarkerFaceColor',frozen_color,'MarkerEdgeColor','k','LineWidth',0.8); % crevasses
% frozen lake frozen
scatter(X_km(1,type_num==5 & Frozen==1),Y_km(1,type_num==5 & Frozen==1),lake_size,...
    's','filled','MarkerFaceColor',frozen_color,'MarkerEdgeColor','k','LineWidth',0.8); % frozen

% surface ice sheet contours
surf_contours = 0:100:1900; % m a.s.l.
[C, h] = contour(BMv5_for_nevis_catchment_noSK.X_m./1e3, BMv5_for_nevis_catchment_noSK.Y_m./1e3, ...
    BMv5_for_nevis_catchment_noSK.S_m, surf_contours,'Color',ice_sheet_contours,'LineWidth',0.8);

% main array GNSS
plot3(xy_sta_22_short(:,1),xy_sta_22_short(:,2),300.*ones(num_sta_22,1),'k^',...
    'MarkerSize',4,'MarkerFaceColor','k','markerEdgeColor','k');

% propagating blister speed radar circles:
for i=1:length(cluster_mid)
plot(squeeze(x_circle(i,2,:)),squeeze(y_circle(i,2,:)),'-','LineWidth',2,'Color',[0.6 0.6 0.6]); 
end

% cluster cOneXiOnS
grey_connections = [0 0 0];
% mid C3
cluster_mid = [168, 173, 183, 185, 201, 223];
cc3 = plot([environs_lakes_2022.X_km(cluster_mid(1)) environs_lakes_2022.X_km(cluster_mid(2))],...
    [environs_lakes_2022.Y_km(cluster_mid(1)) environs_lakes_2022.Y_km(cluster_mid(2))],...
    '-','LineWidth',2,'Color',grey_connections);
cc31 =plot([environs_lakes_2022.X_km(cluster_mid(2)) environs_lakes_2022.X_km(cluster_mid(3))],...
    [environs_lakes_2022.Y_km(cluster_mid(2)) environs_lakes_2022.Y_km(cluster_mid(3))],...
    '-','LineWidth',2,'Color',grey_connections);
cc32 =plot([environs_lakes_2022.X_km(cluster_mid(4)) environs_lakes_2022.X_km(cluster_mid(5))],...
    [environs_lakes_2022.Y_km(cluster_mid(4)) environs_lakes_2022.Y_km(cluster_mid(5))],...
    '-','LineWidth',2,'Color',grey_connections);
cc33 =plot([environs_lakes_2022.X_km(cluster_mid(6)) environs_lakes_2022.X_km(cluster_mid(5))],...
    [environs_lakes_2022.Y_km(cluster_mid(6)) environs_lakes_2022.Y_km(cluster_mid(5))],...
    '-','LineWidth',2,'Color',grey_connections);

lgdB = legend('\sigma_{1}','\sigma_{1} = 200 kPa','24-hr propagation',...
    'HF filling','Moulin filling','Overspill filling',...
    'Crevasses filling','No-exit filling',...
    'HF initial drainage','Moulin initial drainage','Overspill initial drainage','Crevasses initial drainage',...
    'Moulin still draining','Overspill still draining','Crevasses still draining',...
    'HF empty basin','Moulin empty basin','Overspill empty basin','Crevasses empty','No-exit frozen');   
set(lgdB,'Position',[0.885 0.335 0.05 0.04],'FontSize',fontsize-3,'FontName','Avenir',...
    'NumColumns',1,'Box','on','TextColor',lake_blue(1:3))

axis equal; grid on;
date_of_q = sprintf('2022/%d',DOY+1);
text(-44, 15, date_of_q,'FontSize',fontsize+4,'FontName','Helvetica','FontWeight','bold')
set(gca,'FontName','Avenir','FontSize',fontsize,'xtick',-80:10:110,...
    'ytick',-50:10:30,'LineWidth',0.6,'Layer','top')
xlim([plot_x(1) plot_x(2)]); ylim([plot_y(1) plot_y(2)]);

%% print figure
print(gcf,'-dpng','-r300',sprintf('../paperfigs/suppfig2.1_%s_260129.png',cluster_name)); 