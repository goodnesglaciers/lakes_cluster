%% daily lake classifier timeseries movie
% 24-05-01 LAS: make movies of the lake catalogue !
% 24-06-17 LAS: throw out lakes on the nevis boundary, only use lakes that 
% are large enough to HF to bed: call these "HFpossible" lakes; add in
% nevis q model estimates
% 24-06-26 LAS: add in a little timeseries of E over model domain
% 24-07-09 LAS: update with lake-typeing dates based on reviewing all S1 and S2 images
% 24-09-03 LAS: time of drainage plotted as DOY values of 24-hr periods
% 25-04-24 LAS: update FASTER lake locations, two-panel movie w/E+q

clear all; close all;

%% load in needed datasets
% drainage type-ing
% 1 = HF; 2 = moulin; 3 = overspill; 4 = crevasses; 5 = frozen. 6 = slush; 7 = stream; 9 = throw out/inconclusive.
load laketypeing_S1S2MO_2023_251003.mat % collated from Excel sheet
lake_typeing = laketypeing_S1S2_2023;
% lake id = column 1; col 2 = first water in basin; 3 = max lake size;
% 4 = first day of drainage; 5 = still drainin'; 6 = freeze-up date; 
% 7 = type number.
type_num = lake_typeing(:,7); 
lake_num = length(type_num); % number of lakes

% lake locations
load max_extent_indices_2023_GEOD.mat

% bedmap in nevis formation
load('BMv5_for_nevis_catchment_noSK.mat'); % bed machine
H = BMv5_for_nevis_catchment_noSK.S_km - BMv5_for_nevis_catchment_noSK.B_km; % ice thickness
S = BMv5_for_nevis_catchment_noSK.S_km; % ice surface elevation

% station locations
load('polarstereo_stations_2023_short.mat')

% nevis location, ice thickness, and surface/bed for lake locations
%  boundaries and correct dates:
load('environs_lakes_2023B_S1S2MO_251003.mat') %  boundaries and correct dates
environs_lakes_2023B = environs_lakes; 
X_km = environs_lakes_2023B.X_km; Y_km = environs_lakes_2023B.Y_km;

% % within subglacial ROI
% ROI_relevant_xv = [-15 60 104 104 15 -10 -33 -33 -15]; % km in nevis-model space
% ROI_relevant_yv = [15 4 -8 -50 -50 -41 -33 -1 15];
% for i=1:1:length(environs_lakes.H)
%     if inpolygon(environs_lakes.X_km(i),environs_lakes.Y_km(i),...
%             ROI_relevant_xv,ROI_relevant_yv) == 0
%         % throw a NaN
%         environs_lakes.H(i) = NaN;
%         environs_lakes.S(i) = NaN;
%         environs_lakes.X_km(i) = NaN;
%         environs_lakes.Y_km(i) = NaN;
%         environs_lakes.lat(i) = NaN;
%         environs_lakes.lon(i) = NaN;
%         environs_lakes.drainage_type_num(i) = NaN;
%         environs_lakes.laketypeing_dates(i,1:7) = [NaN, NaN, NaN, NaN, NaN, NaN, NaN];
%     else
%     end
% end
% % redefine lake_typeing to just include lakes within nevis domain
% lake_typeing = environs_lakes.laketypeing_dates; % after throwing out lakes on the nevis domain border
% type_num = lake_typeing(:,7); % what type of drainage 
% lake_num = length(type_num); % number of lakes
% environs_lakes_2023 = environs_lakes; 

% load daily_lake_classifier timeseries
load daily_lake_HFpossible_classifier_S1S2MO_2023_251003.mat
daily_lake_classifier = daily_lake_HFpossible_classifier; % name change
daily_lake_classifier(44,189)=1; % still-there HF, on inner S2 image bound
daily_lake_classifier(66,189)=1; % still-there HF, on inner S2 image bound
daily_lake_classifier(103,189)=1; % still-there HF, on inner S2 image bound
daily_lake_classifier(131,189)=1; % still-there HF, on inner S2 image bound
% % within subglacial ROI:
% for i=1:1:length(environs_lakes_2023.S)
%     if isnan(environs_lakes_2023.S(i))
%         % throw a NaN
%         daily_lake_classifier(i,:) = NaN;
%     else
%     end
% end

% load nevis discharge q outputs
    addpath('nevis')
    % load initial timestep
    fn = 'nevis_lakesix_2023_noSK_300m_ub';
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
        
    %boundary curve
    if ~isempty(gg.n1)
    x_out = gg.nx(gg.n1); y_out = gg.ny(gg.n1);
    tmp = nevis_orderboundary(x_out,y_out); x_out = x_out(tmp); y_out = y_out(tmp); % reorder to follow boundary
    else x_out = []; y_out = [];
    end

%% FIGURES
S=25; % marker size

%% Lake-drainage Type. 2023 Sentinel-2. no change in time. 
figure(1); clf;
Fig1 = set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',1.0.*[1 1 33.3 17.3]);
axe1 = axes('Position',[0.05 0.12 0.90 0.80],'Box','on','NextPlot','add');
hold on; 
scatter(X_km(1,type_num==1), Y_km(1,type_num==1),'dr','filled','MarkerEdgeColor','k'); % HF
scatter(X_km(1,type_num==2), Y_km(1,type_num==2),'^b','filled','MarkerEdgeColor','k'); % moulin
scatter(X_km(1,type_num==3), Y_km(1,type_num==3),'o','filled','MarkerFaceColor',[0.804 0.704 0.804]); % overspill
scatter(X_km(1,type_num==4), Y_km(1,type_num==4),'vy','filled','MarkerEdgeColor','k'); % crevasses
plot(X_km(1,type_num==5), Y_km(1,type_num==5),'Xk','markerSize',7); % frozen
title('Daily Lake Classification. 2023 Sentinel-2. HF possible.')
ylabel('North–South [ km ]')
xlabel('East–West [ km ]')
lgd=legend('hydro-fracture','moulin','overspill','water-filled crevasses','no-exit frozen');
set(lgd,'NumColumns',1,'Location','NorthWest')
axis equal 
%ylim([0.1e4 2.4e4]); xlim([0 4e4]) % pixels
% ylim([7565 7665]); xlim([520 660]) % UTM 22N
ylim([-55 30]); xlim([-60 120]) % nevis model domain
box on
grid on
set(gca,'LineWidth',0.8,'FontName','Avenir')
%print(gcf,'-dpng','-r300','classifier_2023_HFpossible_type_240830.png'); 

%% load in nevis outputs for Q_out, basal melt, and E across domain
for j=1
    disp(['Frame ',num2str(j),' / ',num2str(length(j)),' ...']);
    % load timestep
    load([fn,'/',int2four(j)]);
    % extract new variables
    if isfield(vv,'nbdy'), gg = nevis_label(gg,vv.nbdy); aa.phi = aa.phi_a(gg.nbdy); end % remask
    aa = nevis_inputs(vv.t,aa,pp,gg,oo);
    oo.evaluate_variables = 1; oo.evaluate_residual = 0; oo.evaluate_jacobian = 0; [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
    vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
    nevis_unpack(aa,vv2);

    % append to model-domain Q timeseries
    time = (ps.t/(24*60*60))*[tt.t];
    Q_in = ps.Q*[tt.Q_in];
    Q_out = ps.Q*[tt.Q_out];
    basal_melt = (ps.m*ps.x^2)*[tt.m];
    surface_runoff = (ps.m*ps.x^2)*[tt.E];

    % clear this timestep
    clear tt vv vv2
end
% 
% for j=2:300
%     disp(['Frame ',num2str(j),' / ',num2str(length(j)),' ...']);
%     % load timestep
%     load([fn,'/',int2four(j)]);
%     % extract new variables
% %     if isfield(vv,'nbdy'), gg = nevis_label(gg,vv.nbdy); aa.phi = aa.phi_a(gg.nbdy); end % remask
% %     aa = nevis_inputs(vv.t,aa,pp,gg,oo);
% %     oo.evaluate_variables = 1; oo.evaluate_residual = 0; oo.evaluate_jacobian = 0; [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
% %     vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
% %     nevis_unpack(aa,vv2);
%     
%     % append to model-domain Q timeseries
%     time = horzcat(time,(ps.t/(24*60*60))*[tt.t]);
%     Q_in = horzcat(Q_in,ps.Q*[tt.Q_in]);
%     Q_out = horzcat(Q_out,ps.Q*[tt.Q_out]);
%     basal_melt = horzcat(basal_melt,(ps.m*ps.x^2)*[tt.m]);
%     surface_runoff = horzcat(surface_runoff,(ps.m*ps.x^2)*[tt.E]);
% 
%     % clear this timestep
%     clear tt %vv vv2
% end
% 
% % find which rows are unique in tt
% [unique_time, ia] = unique(time);
% % Keep unique rows
% time_keep = time(ia);
% Q_in_keep = Q_in(ia);
% Q_out_keep = Q_out(ia);
% basal_melt_keep = basal_melt(ia);
% surface_runoff_keep = surface_runoff(ia);
% % store and save
% nevis_lakesix_2023_noSK_300m_ub_unique.time = time_keep;
% nevis_lakesix_2023_noSK_300m_ub_unique.Q_in = Q_in_keep;
% nevis_lakesix_2023_noSK_300m_ub_unique.Q_out = Q_out_keep;
% nevis_lakesix_2023_noSK_300m_ub_unique.basal_melt = basal_melt_keep;
% nevis_lakesix_2023_noSK_300m_ub_unique.surface_runoff = surface_runoff_keep;
% save nevis_lakesix_2023_noSK_300m_ub_unique.mat nevis_lakesix_2023_noSK_300m_ub_unique
load nevis_lakesix_2023_noSK_300m_ub_unique.mat

% axes(axe9)
% 
%     % N saved at lake and GPS locations
%     colormap1=bone(4);
%     colormapt=parula(22);
%     xlim([140 330]); ylim([-1 3]);
%     hold on;
%     %text(226, 2.5, 'd.', 'FontSize',m+1,'FontWeight','bold');
%     ylabel('N [ MPa ]'); hold on;
%     xlabel('DOY 2022 (Model)');
%     test=(ps.phi/10^6).*[tt.pts_N];
%     for i=1:2:3; plot([tt.t].*10,test(i,:),'-','Color',colormap1(i,:),'LineWidth',1.1); end
%     for i=5:3:22; plot([tt.t].*10,test(i,:),'-','Color',colormapt(i,:),'LineWidth',1.1); end
%     rectangle("Position",[lake_drain_days_2022(1) -1 1 4],'EdgeColor','none','FaceColor',[1 0.85 0.85])
%     rectangle("Position",[lake_drain_days_2022(2) -1 1 4],'EdgeColor','none','FaceColor',[1 0.85 0.85])
%     rectangle("Position",[lake_drain_days_2022(3) -1 1 4],'EdgeColor','none','FaceColor',[1 0.85 0.85])
%     plot([doy doy], [-5  5], '-k', 'LineWidth', 1);
%     plot([1 365], [0 0], '--k', 'LineWidth', 0.8);
%     for i=1:2:3; plot([tt.t].*10,test(i,:),'-','Color',colormap1(i,:),'LineWidth',1.1); end
%     for i=5:3:22; plot([tt.t].*10,test(i,:),'-','Color',colormapt(i,:),'LineWidth',1.1); end
%     set(gca, 'FontSize',mm)
%     load('../RACMO/station_names.mat')
%     lgdnames=[station_names(1:2:3)';station_names(5:3:22)'];
%     %legend(lgdnames,'Location','EastOutside');
%     celllgd = cellstr(lgdnames);
%     %columnlegend(3, celllgd,'Location','SouthEast'); %,'Position',[-1 232 1 1]);
%     %columnlegend(3, celllgd,'Position',[0.8 -10000 10 10]);
%    
%     legend(celllgd, 'Location','NorthEast','NumColumns',[4]);
%     legend boxoff
%     
%     set(gca, 'FontSize',mm)
%     set(gca,'TickDir','out','LineWidth',1.01); box on;


%% daily catalogue -- q base
time_days = 1:1:365;
TriangleSize = 5; 
HF_size = 70;
lake_size = 30;
mm=10;
% load cmapland2.mat

filling_color = [0.2 0.6 1.0];
draining_color = [1.0 0.8 0.0];
still_draining_color = [0.2 0.8 0.2];
frozen_color = [0.8 0.8 0.8];

for j=150:1:300 %length(time_days)

Fig2 = figure(j); clf;
set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',1.0.*[1 1 33.3 17.3*2]);
axe1 = axes('Position',[0.05 0.58 0.90 0.43],'Box','on','NextPlot','add','XTickLabel',[]);
axe2 = axes('Position',[0.05 0.20 0.90 0.43],'Box','on','NextPlot','add');
axe3 = axes('Position',[0.05 0.035 0.45 0.15],'Box','on','NextPlot','add');

axes(axe1); 
hold on; 

% nevis runoff: plot frames
    disp(['Frame ',num2str(j),' / ',num2str(length(j)),' ...']);
    % load timestep
    load([fn,'/',int2four(j)]);
    % extract new variables
    if isfield(vv,'nbdy'), gg = nevis_label(gg,vv.nbdy); aa.phi = aa.phi_a(gg.nbdy); end % remask
    aa = nevis_inputs(vv.t,aa,pp,gg,oo);
    oo.evaluate_variables = 1; oo.evaluate_residual = 0; oo.evaluate_jacobian = 0; [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
    vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
    nevis_unpack(aa,vv2);
    clear tt vv vv2

    % runoff plotting
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny;
    zz = (ps.m*10^3*pd.td)*reshape(E,gg.nI,gg.nJ); 
    cax = [0 80]; 
    caxis(cax); 
    axis image;
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis
    load cmapq2; cmap2=cmap; colormap(gca, cmap2); %cmap2(64,1:3) = NaN;
    hand = pcolor(xx,yy,zz); % shading interp;
    set(hand,'linestyle','none','FaceAlpha',1.0); 
    % color bar
    cbar_size = [0.75 0.17 0.20 0.01];
    cx5 = colorbar('peer',gca,'horizontal','position',cbar_size,'XAxislocation','bottom');
    tmp = get(cx5,'XTick');  labels = tmp;  
    xlabel(cx5, 'Surface Runoff [ mm d^{-1} ]')
    set(cx5,'XTickLabel',labels,'tickdir','out');
    
% surface elevation contours
[surf_contours,surf_labels] = contour(BMv5_for_nevis_catchment_noSK.X_km, BMv5_for_nevis_catchment_noSK.Y_km, ...
    double(BMv5_for_nevis_catchment_noSK.S_km),0:100:1900,'LineWidth',0.8,...
    'LineStyle','-','EdgeColor',[0.7 0.7 0.7]);
clabel(surf_contours,surf_labels,'FontName','Avenir','Color',[0.4 0.4 0.4],...
    'FontSize',mm-2,'labelspacing', 700)

% gps stations
plot(xy_sta_23_short(:,1),xy_sta_23_short(:,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');

% lake filling
Filling = daily_lake_classifier(:,j)==1;
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
Draining = daily_lake_classifier(:,j)==2;
% HF lake draining
scatter(X_km(1,type_num==1 & Draining==1),Y_km(1,type_num==1 & Draining==1),HF_size+30,...
    'p','filled','MarkerFaceColor',draining_color,'MarkerEdgeColor','k'); % HF
% moulin lake draining
scatter(X_km(1,type_num==2 & Draining==1),Y_km(1,type_num==2 & Draining==1),lake_size+30,...
    'v','filled','MarkerFaceColor',draining_color,'MarkerEdgeColor','k'); % moulin
% overspill lake draining
scatter(X_km(1,type_num==3 & Draining==1),Y_km(1,type_num==3 & Draining==1),lake_size+30,...
    'o','filled','MarkerFaceColor',draining_color,'MarkerEdgeColor','k'); % overspill
% crevasses draining
scatter(X_km(1,type_num==4 & Draining==1),Y_km(1,type_num==4 & Draining==1),lake_size+30,...
    'd','filled','MarkerFaceColor',draining_color,'MarkerEdgeColor','k'); % crevasses

% lake draining -- continued
Still_Draining = daily_lake_classifier(:,j)==3;
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
Empty = daily_lake_classifier(:,j)==3;
% HF empty basin
scatter(X_km(1,type_num==1 & Empty==1),Y_km(1,type_num==1 & Empty==1),HF_size,...
    'p','filled','MarkerFaceColor',frozen_color,'MarkerEdgeColor','k','LineWidth',0.8); % HF

% lake frozen
Frozen = daily_lake_classifier(:,j)==4;
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

title(sprintf('2023/%d',time_days(j)), "FontName",'AvenirBold','FontWeight','bold',FontSize=mm+1);
ylabel('North–South [ km ]'); %xlabel('East–West [ km ]')

lgd=legend('Surface Runoff {\itE}','Surface Elevation [ m a.s.l. ]','GNSS station',...
    'Hydro-fracture (HF) filling','Moulin (MU) filling','Overspill (OS) filling','Crevasses (CV) filling','No-exit filling', ...
    'HF draining','MU initial drainage','OS initial drainage','CV initial drainage',...
    'MU still draining','OS still draining','CV still draining',...
    'HF empty basin',...
    'Moulin dry/frozen','Overspill dry/frozen','Crevasses dry/frozen','No-exit frozen');
set(lgd,'NumColumns',1,'Location','NorthWest','FontSize',mm-2)
axis equal 
ylim([-50.05 28.25]); xlim([-80 104]) % nevis model domain
box on; grid on
set(gca,'LineWidth',1.0,'FontName','Avenir','FontSize',mm)
set(gca,'layer','top')

axes(axe2); 
hold on; 

% nevis discharge: plot frames
    disp(['Frame ',num2str(j),' / ',num2str(length(j)),' ...']);
    % load timestep
    load([fn,'/',int2four(j)]);
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
    load cmapq2; cmap2=cmap; colormap(gca, cmap2); %cmap2(64,1:3) = NaN;
    hand = pcolor(xx,yy,log10(zz)); % shading interp;
    set(hand,'linestyle','none','FaceAlpha',1.0); 
    caxis(log10(cax)); 
    axis image;
    % color bar
    cbar_size = [0.75 0.10 0.20 0.01];
    cx5 = colorbar('peer',gca,'horizontal','position',cbar_size,'XAxislocation','bottom');
    tmp = get(cx5,'XTick');  labels = 10.^tmp;  
    xlabel(cx5, 'Subglacial Discharge {\itq} [ m^2 s^{-1} ]')
    set(cx5,'XTickLabel',labels,'tickdir','out');
    
% surface elevation contours
[surf_contours,surf_labels] = contour(BMv5_for_nevis_catchment_noSK.X_km, BMv5_for_nevis_catchment_noSK.Y_km, ...
    double(BMv5_for_nevis_catchment_noSK.S_km),0:100:1900,'LineWidth',0.8,...
    'LineStyle','-','EdgeColor',[0.7 0.7 0.7]);
clabel(surf_contours,surf_labels,'FontName','Avenir','Color',[0.4 0.4 0.4],...
    'FontSize',mm-2,'labelspacing', 700)

% gps stations
plot(xy_sta_23_short(:,1),xy_sta_23_short(:,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');

% lake filling
Filling = daily_lake_classifier(:,j)==1;
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
Draining = daily_lake_classifier(:,j)==2;
% HF lake draining
scatter(X_km(1,type_num==1 & Draining==1),Y_km(1,type_num==1 & Draining==1),HF_size+30,...
    'p','filled','MarkerFaceColor',draining_color,'MarkerEdgeColor','k'); % HF
% moulin lake draining
scatter(X_km(1,type_num==2 & Draining==1),Y_km(1,type_num==2 & Draining==1),lake_size+30,...
    'v','filled','MarkerFaceColor',draining_color,'MarkerEdgeColor','k'); % moulin
% overspill lake draining
scatter(X_km(1,type_num==3 & Draining==1),Y_km(1,type_num==3 & Draining==1),lake_size+30,...
    'o','filled','MarkerFaceColor',draining_color,'MarkerEdgeColor','k'); % overspill
% crevasses draining
scatter(X_km(1,type_num==4 & Draining==1),Y_km(1,type_num==4 & Draining==1),lake_size+30,...
    'd','filled','MarkerFaceColor',draining_color,'MarkerEdgeColor','k'); % crevasses

% lake draining -- continued
Still_Draining = daily_lake_classifier(:,j)==3;
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
Empty = daily_lake_classifier(:,j)==3;
% HF empty basin
scatter(X_km(1,type_num==1 & Empty==1),Y_km(1,type_num==1 & Empty==1),HF_size,...
    'p','filled','MarkerFaceColor',frozen_color,'MarkerEdgeColor','k','LineWidth',0.8); % HF

% lake frozen
Frozen = daily_lake_classifier(:,j)==4;
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

%title(sprintf('2023/%d',time_days(j)), "FontName",'Avenir','FontWeight','bold',FontSize=mm+1);
ylabel('North–South [ km ]'); xlabel('East–West [ km ]')

lgd=legend('Subglacial discharge {\itq}');
set(lgd,'NumColumns',1,'Location','NorthWest','FontSize',mm-2)
axis equal 
ylim([-50.05 28.25]); xlim([-80 104]) % nevis model domain
box on; grid on
set(gca,'LineWidth',1.0,'FontName','Avenir','FontSize',mm)
set(gca,'layer','top')

% timeseries of E integrated across the domain
axes(axe3)
hold on; 
% model output
    plot(nevis_lakesix_2023_noSK_300m_ub_unique.time, nevis_lakesix_2023_noSK_300m_ub_unique.surface_runoff, '-b','LineWidth', 1.3);
    plot(nevis_lakesix_2023_noSK_300m_ub_unique.time, nevis_lakesix_2023_noSK_300m_ub_unique.basal_melt, ':b','LineWidth', 1.3);
    plot(nevis_lakesix_2023_noSK_300m_ub_unique.time, nevis_lakesix_2023_noSK_300m_ub_unique.Q_out, '-','LineWidth', 1.4,'Color',draining_color);
% moving time marker
plot([j j],[1 7999],'-','LineWidth',1.5,'Color',[0.4940 0.1840 0.5560]);
plot(j,6000,'v','MarkerSize',TriangleSize+2,'MarkerFaceColor',[0.4940 0.1840 0.5560],'Color',[0.4940 0.1840 0.5560]);
% model output
    plot(nevis_lakesix_2023_noSK_300m_ub_unique.time, nevis_lakesix_2023_noSK_300m_ub_unique.Q_out, '-','LineWidth', 1.4,'Color',draining_color);
    plot(nevis_lakesix_2023_noSK_300m_ub_unique.time, nevis_lakesix_2023_noSK_300m_ub_unique.surface_runoff, '-b','LineWidth', 1.3);

    ylabel('Q [ m^{3} s^{-1} ]'); xlabel('Day of Year, 2023');
    
    lgdnames2=[{'Surface Runoff'};{'Basal Melt'};{'Proglacial Discharge'}];
    celllgd2 = cellstr(lgdnames2);
    

xlim([150 300]); ylim([0 6000]); 
set(gca,'xtick',[150:25:300],'LineWidth',1.0,'FontName','Avenir','FontSize',mm)
box on; grid on

lgd3 = columnlegend(1, celllgd2,'Location','NorthEast','FontSize',mm-2); 

%% mini Green Greenland
lake_blue = [0    0.4470    0.7410   0.08]; % fourth value is FaceAlpha built in!
goldenrod = [1.0 0.8 0.0];
axe_map = axes('Position',[0.005 0.39 0.19 0.19],'Box','off','XTickLabel',[],'YTickLabel',[]);
axes(axe_map)
ax = worldmap('greenland');
greenland = shaperead('landareas', 'UseGeoCoords', true,...
  'Selector',{@(name) strcmp(name,'Greenland'), 'Name'});
patchm(greenland.Lat,greenland.Lon,'FaceColor', [0.4 0.6 0.4],'EdgeColor',lake_blue(1:3),...
    'LineWidth',1.1); hold on;
%plotm([69],[-49],'p','MarkerSize',20,'MarkerFaceColor',lake_blue(1:3),'MarkerEdgeColor',[1 1 1]); 
fillm([68.2 69.2 69.2 68.2],[-47.5 -47.5 -50.5 -50.5],'w','LineWidth',1.1,'FaceColor',goldenrod,'EdgeColor',goldenrod)
axis off; framem off; gridm off; mlabel off; plabel off; 

% print figure
%print(gcf,'-dpng','-r300',sprintf('classifier_HFpossible_E_q_base/classifier_2023_HFpossible_E_q_Q_250424_%d.png',time_days(j))); 
%close all

end


%% ffmpeg line
%% ffmpeg -r 6 -f image2 -pattern_type glob -i "classifier*.png" -vcodec libx264 -crf 25 -vf "crop=trunc(iw/2)*2:trunc(ih/2)*2" -pix_fmt yuv420p output-2023-classifier-bedtopo.mp4
%  Friendly for large image sizes:
 % ffmpeg -r 6 -f image2 -pattern_type glob -i "classifier*.png" \
 %  -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" \
 %  -vcodec libx264 -crf 20 -pix_fmt yuv420p \
 %  -profile:v high -level 4.1 -x264opts keyint=6:min-keyint=6:no-scenecut \
 %  output-qt-friendly.mp4
