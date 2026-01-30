%% clusters - homogeneous Poisson process 
% 24-05-06 LAS -- HF drainages in 2022 from S2
% 24-05-14 LAS -- automatically identify clusters
% 24-05-15 LAS -- HF, MU, OS all together in 2022 from S2
% 24-06-17 LAS -- "HFpossible" candidate lakes
% 24-06-24 LAS -- order lakes by their elevation
% 24-08-23 LAS -- slim panels down for paper figures :D -- just HF
% 24-08-26 LAS -- dates now using S1+S2 images in 2022 and 2023
% 24-09-03 LAS -- time of drainage plotted as DOY values of 24-hr periods
% 24-09-03 LAS -- Change question to ask what's stat sig in the drainage
% catchment? (i.e., just the region of the ROI that's relevant to the GNSS array)
% 24-09-06 LAS -- racmo at HF event
% 24-10-23 LAS -- update cut-off lakes in C1 and C4 cluster on edge of S2
% (L34) image
% 25-02-13 LAS -- drainage-mechanisms by elevation bar plots, new panels
% 25-03-16 LAS -- add a location map!
% 25-04-16 LAS -- update with FASTER lake boundaries, locations, elevations
% 25-05-28 LAS -- v_crit and lake formation (in DOY)
% 26-01-28 LAS -- repository checks

clear all; close all;

%% CLASSIFIER
% 1 = HF; 2 = moulin; 3 = overspill; 4 = crevasses; 5 = frozen. 9 = throw out/inconclusive.
load('daily_lake_HFpossible_classifier_S1S2_2022_240829.mat'); % S1+S2 images
daily_lake_classifier(34,196) = 3; % drained HF, L34 on inner S2 image bound

% lakes large enough to HF and in nevis domain:
% load lake boundaries, centre points, drainage dates from FASTER for 2022
load('environs_lakes_2022B_250416.mat') % drainage dates & boundaries
environs_lakes_2022 = environs_lakes; 
lake_typeing = environs_lakes.laketypeing_dates; % has NaN outside the nevis area
lake_typeing(environs_lakes.lake_HF_possible_days==0,7) = NaN; % 0 days of HF possible (i.e., lake is too small to HF)
lake_elev = environs_lakes.S; % lake elevation from BedMachine

% load available image dates for FASTER for 2022
load('all_JD_nums_2022.mat'); % all julian day numbers of S2 images
load('imageID_2022_250523.mat'); % JD numbers of S2 images available for each specific lake

% ROI that's relevant to the GNSS array 
% if the lake location is not within relevant catchment, throw it out !!! 
% after FASTER locations update (2025-April):
ROI_relevant_xv = [-15 60 104 104 15 -10 -33 -33 -15]; % X km in nevis-model space
ROI_relevant_yv = [15 4 -8 -50 -50 -41 -33 -1 15]; % Y km in nevis-model space
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

% lake_typeing of lakes within the ROI (for elevation-bin analysis 25-02-13)
lake_typeing_all = lake_typeing;
% count up lakes by elevation bands
elev_band_interval = 100;
elev_bands = 700:elev_band_interval:1900; % 100 m a.s.l.
for i=1:1:length(elev_bands)-1
    % find lakes within this elevation band
    counting_elev_lakes(i).id = find(lake_elev>elev_bands(i) & lake_elev<=elev_bands(i+1));
    % count how many HF/moulin/overspill/frozen events there are
    counting_elev_lakes(i).mechanism(1,1) = length(find(lake_typeing_all(counting_elev_lakes(i).id,7)==1)); % HF
    counting_elev_lakes(i).mechanism(1,2) = length(find(lake_typeing_all(counting_elev_lakes(i).id,7)==2)); % in-lake moulin
    counting_elev_lakes(i).mechanism(1,3) = length(find(lake_typeing_all(counting_elev_lakes(i).id,7)==3)); % overspill
    counting_elev_lakes(i).mechanism(1,4) = length(find(lake_typeing_all(counting_elev_lakes(i).id,7)==5)); % frozen
    % count how many lakes there are
    counting_elev_lakes(i).number_form = sum(counting_elev_lakes(i).mechanism(1,1:4));
    % percentage of HF out of all lakes that form
    counting_elev_lakes(i).percent_form = 100.*(counting_elev_lakes(i).mechanism./counting_elev_lakes(i).number_form);
    % count how many lakes drain
    counting_elev_lakes(i).number_drain = sum(counting_elev_lakes(i).mechanism(1,1:3));
    % percentage of HF out of drainage events (i.e., lakes that drain)
    counting_elev_lakes(i).percent_drain = 100.*(counting_elev_lakes(i).mechanism(1:3)./counting_elev_lakes(i).number_drain);
    % save elevation band lower bound
    counting_elev_lakes(i).elev_lower_bound = elev_bands(i);
end
counting_elev_lakes_ROI_22 = counting_elev_lakes;
save counting_elev_lakes_ROI_22.mat counting_elev_lakes_ROI_22

% % polygon ROI check
% figure; plot(ROI_relevant_xv,ROI_relevant_yv)
% hold on; plot(environs_lakes.X_km,environs_lakes.Y_km,'o')

% RACMO runoff at time of HF
load cumulative_runoff_last_day_to_HF_2022.mat

% preliminaries for homo poisson
alpha = 0.05; % significance threshold (one tailed)
t = 1:1:70; % time window [ days ]
m = 0:1:250; % number of events [ cluster size ]

%% drainage dates for each drainage type
% sort for HF drainages

% HF
id_HF_drainages = find(lake_typeing(:,7)==1); % HF drainages row IDs [ HF == 1 ]
save id_HF_drainages_2022_250415.mat id_HF_drainages
days_HF_drainages = daily_lake_classifier(id_HF_drainages,:); % HF drainages timeseries
days_HF_drainages(days_HF_drainages~=2) = NaN; % NaN if not in drainage window
days_HF_drainages(days_HF_drainages==2) = 1; % drainage window day ---> 1

% % get the in-between days, and not the initial day, --> 1
days_HF_drainages2 = NaN(length(id_HF_drainages),365); % prefill with NaNs
for i=1:1:length(id_HF_drainages)
    % if just one day draining window: it's that day
    if environs_lakes.laketypeing_dates(id_HF_drainages(i),3) == environs_lakes.laketypeing_dates(id_HF_drainages(i),4)
        days_HF_drainages2(i,environs_lakes.laketypeing_dates(id_HF_drainages(i),4)) = 1;
    elseif environs_lakes.laketypeing_dates(id_HF_drainages(i),3)+1 == environs_lakes.laketypeing_dates(id_HF_drainages(i),4)
        days_HF_drainages2(i,environs_lakes.laketypeing_dates(id_HF_drainages(i),4)) = 1;
    else 
        days_HF_drainages2(i,(environs_lakes.laketypeing_dates(id_HF_drainages(i),3)+1):environs_lakes.laketypeing_dates(id_HF_drainages(i),4)) = 1;
   end
end
%figure(1); imagesc(days_HF_drainages); colorbar; xlim([160 240]); title('HF start stop'); shg
%figure(2); imagesc(days_HF_drainages2); colorbar; xlim([160 240]); title('HF active drainage 24-hr'); shg
% calculate number of HF lakes draining on all days
days_HF_drainages_sum = nansum(days_HF_drainages2,1);
% calculate error in each lake-drainage date (this is the window size)
days_HF_drainages_window = nansum(days_HF_drainages2,2);
average_HF_window = nanmean(days_HF_drainages_window);
% drainages
drainages = length(id_HF_drainages) % number of HF drainages in 2022
drainages_start = 185; % [ DOY ]
drainages_stop = 215; % [ DOY ]
days = drainages_stop-drainages_start; % 24-hr windows HF drainages occur in in 2022
lambda = drainages./days; % [ events / day ]
% lake elevations
elevations_HF = lake_elev(id_HF_drainages);

%% %%%%%%%%% v_crit and lake formation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get v_crit and lake-formation time windows
days_HF_form = daily_lake_classifier(id_HF_drainages,:); % HF drainages timeseries
for i=1:length(days_HF_form(:,1))
    % what's the first image we have of this lake?
    day_lake_form_last(i,1) = lake_typeing(id_HF_drainages(i,1),2); % image date of lake first formed
    index1 = find(all_JD_nums(:,1)==day_lake_form_last(i,1)); % Sentinel image index on that day
    index2 = find(imageID_2022.lake(id_HF_drainages(i),1).plot_IDs==index1); % Sentinel image on day lake first formed (lake specific)
    % what's the last image available before we see a lake?  (date of last dry basin) 
    index3 = imageID_2022.lake(id_HF_drainages(i),1).plot_IDs(index2-1); % image index before lake first formed (lake specific)
    day_lake_form_first(i,1) = all_JD_nums(index3,1)+1; % date one day after image of dry lake basin (lake specific)
    error_day_lake_form(i,1) = day_lake_form_last(i,1)-day_lake_form_first(i,1); % days of spread (count 24-hr periods)
    % v_crit model estimate
    day_HF_possible(i,1) = find(days_HF_form(i,:)==1,1,'first'); % first day of lake being HFpossible
end
% for plotting
formation_average = (day_lake_form_first+day_lake_form_last)/2; % midpoint [ DOY ]
formation_error = day_lake_form_last-formation_average; % error [ days ]

%% HF in 2022 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf; 
set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',2.1.*[0 0 8.3 13]);
axe1 = axes('Position',[0.085 0.54 0.70 0.45],'Box','on','NextPlot','add');
axe2 = axes('Position',[0.095 0.815 0.25 0.17],'Box','on','NextPlot','add');
axe3 = axes('Position',[0.790 0.54 0.18 0.45],'Box','on','NextPlot','add');

% plotting preliminaries
fontsize = 9; 
TriangleSize = 7;
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

% hydro-fracture catalogue 2022 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(axe1); hold on 
yyaxis left

i=2; % ordering is for populating legend
errorbar(formation_average, elevations_HF,[],[],formation_error,formation_error,'k.',...
    'LineWidth',0.7,'Marker','none','CapSize',3); % lake formation
plot(day_HF_possible,elevations_HF,'dk',...
    'markerFaceColor','w','MarkerSize',5); % v_crit
plot(1:1:365,days_HF_drainages2(i,:)+elevations_HF(i)-days_HF_drainages2(i,:),...
        'dk','markerFaceColor',goldenrod,'MarkerSize',5) % HF

% RACMO : passing melt thresholds at elevation bands
load('runoff_elev_2022_300m.mat')
melt_threshold = [50, 250, 450, 650, 850, 1050, 1250]; % [mm]
for i=1:1:length(runoff_elev_2022_300m.elev_bands)-1
    for j=1:length(melt_threshold)
    if runoff_elev_2022_300m.cumulative(i,end) >= melt_threshold(j)
    DOY_threshold(i,j) = find(runoff_elev_2022_300m.cumulative(i,:)>melt_threshold(j),1,'first'); % first day above 0.25 m
    else
    DOY_threshold(i,j) = 1e9;
    end
   end
end
melt_colors = [[0.9 0.9 0.9],;
       [0.85 0.85 0.85],;
       [0.8 0.8 0.8],;
       [0.75 0.75 0.75],;
       [0.7 0.7 0.7],;
       [0.65 0.65 0.65],;
       [0.5 0.5 0.5]];
area(DOY_threshold(:,1),runoff_elev_2022_300m.elev_bands(1:end-1),'EdgeColor','none','FaceColor',melt_colors(1,:)); 
area(DOY_threshold(:,2),runoff_elev_2022_300m.elev_bands(1:end-1),'EdgeColor','none','FaceColor',melt_colors(2,:)); 
area(DOY_threshold(:,3),runoff_elev_2022_300m.elev_bands(1:end-1),'EdgeColor','none','FaceColor',melt_colors(3,:)); 
area(DOY_threshold(:,4),runoff_elev_2022_300m.elev_bands(1:end-1),'EdgeColor','none','FaceColor',melt_colors(4,:)); 
area(DOY_threshold(:,5),runoff_elev_2022_300m.elev_bands(1:end-1),'EdgeColor','none','FaceColor',melt_colors(5,:)); 
area(DOY_threshold(:,6),runoff_elev_2022_300m.elev_bands(1:end-1),'EdgeColor','none','FaceColor',melt_colors(6,:));
area(DOY_threshold(:,7),runoff_elev_2022_300m.elev_bands(1:end-1),'EdgeColor','none','FaceColor',melt_colors(7,:));

% HF possible vcrit
errorbar(formation_average, elevations_HF,[],[],formation_error,formation_error,'k.',...
    'LineWidth',0.7,'Marker','none','CapSize',3); % lake formation
plot(day_HF_possible,elevations_HF,'dk',...
    'markerFaceColor','w','MarkerSize',5); % v_crit

% HF clusters
rectangle('Position',[213.5 1100 2 135],'EdgeColor',lake_blue(1:3),'FaceColor','none','LineWidth',1.1);
rectangle('Position',[206.5 1155 4 150],'EdgeColor',lake_blue(1:3),'FaceColor','none','LineWidth',1.1);
rectangle('Position',[194.5 810 1 160],'EdgeColor',lake_blue(1:3),'FaceColor','none','LineWidth',1.1);
for i=1:1:length(days_HF_drainages2(:,1))
    plot(1:1:365,days_HF_drainages2(i,:)+elevations_HF(i)-days_HF_drainages2(i,:),...
        'k-','LineWidth',0.9) 
end
for i=1:1:length(days_HF_drainages2(:,1))
    plot(1:1:365,days_HF_drainages2(i,:)+elevations_HF(i)-days_HF_drainages2(i,:),...
        'dk','markerFaceColor',goldenrod,'MarkerSize',5) 
end
xlim([175 216]); ylim([700 1800]);
set(gca,'ytick',500:100:1800);
set(gca,'LineWidth',1.1,'tickdir','in','TickLength',[0 0],'FontSize',fontsize,'FontName','Avenir');
set(gca,'ycolor','k');
ylabel('Ice-sheet Surface Elevation [ m a.s.l. ]','FontSize',fontsize);
xlabel('Day of Year, 2022','FontSize',fontsize)
set(gca, 'SortMethod', 'depth')
text(148.5, 1800, 'a','FontName','Arial','FontWeight','bold','FontSize',fontsize+4);
% % instrumented lakes (text too busy for paperfig)
% text(192.25, elevations_HF(9),'L1A','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10); % 81
% text(195.75, elevations_HF(8), 'L1B','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10); % 72
% text(211.3+0.01, elevations_HF(15), 'L2A','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10);
% text(211.3+0.01, elevations_HF(12), 'L2B','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10);
% text(211.3+0.01, elevations_HF(20), 'L2C','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10);
% text(211.3+0.01, elevations_HF(14), 'L2D','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10);
% % % other lakes
% text(206.5+1, elevations_HF(16)-30, 'L168','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10);
% text(find(days_HF_drainages2(1,:)==1,1,'first'), elevations_HF(1), 'L9','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages2(2,:)==1,1,'first')+5.25, elevations_HF(2), 'L21','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% %text(find(days_HF_drainages2(3,:)==1,1,'first'), elevations_HF(3), 'L26','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages2(3,:)==1,1,'first')-0.2, elevations_HF(3), 'L34','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages2(4,:)==1,1,'first'), elevations_HF(4), 'L39','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages2(5,:)==1,1,'first')+3.25, elevations_HF(5), 'L47','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages2(6,:)==1,1,'first'), elevations_HF(6), 'L63','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages2(7,:)==1,1,'first')-0.2, elevations_HF(7), 'L72','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% % text(find(days_HF_drainages2(11,:)==1,1,'first'), elevations_HF(11), 'L85','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% % text(find(days_HF_drainages2(12,:)==1,1,'first'), elevations_HF(12), 'L86','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages2(10,:)==1,1,'first'), elevations_HF(10), 'L102','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages2(11,:)==1,1,'first'), elevations_HF(11), 'L108','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% % text(find(days_HF_drainages2(15,:)==1,1,'first'), elevations_HF(15), 'L118','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% % text(find(days_HF_drainages2(16,:)==1,1,'first'), elevations_HF(16), 'L137','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages2(13,:)==1,1,'first'), elevations_HF(13), 'L146','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages2(17,:)==1,1,'first')-0.2, elevations_HF(17)-30, 'L173','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages2(18,:)==1,1,'first'), elevations_HF(18), 'L183','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% % text(find(days_HF_drainages2(24,:)==1,1,'first'), elevations_HF(24), 'L184','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages2(19,:)==1,1,'first')-0.2, elevations_HF(19), 'L185','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages2(21,:)==1,1,'first'), elevations_HF(21)-5, 'L201','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% % text(find(days_HF_drainages2(28,:)==1,1,'first'), elevations_HF(28)-15, 'L209','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% % text(find(days_HF_drainages2(29,:)==1,1,'first'), elevations_HF(29)+5, 'L214','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% % text(find(days_HF_drainages2(30,:)==1,1,'first'), elevations_HF(30)+10, 'L220','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages2(22,:)==1,1,'first'), elevations_HF(22)+10, 'L223','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% cluster labels
grey_connections = lake_blue(1:3).*1.3; %[0.7 0.7 0.7];
text(194.5, 1005,'C1','FontSize',fontsize+6,'FontName','Avenir','FontAngle','italic',...
    'FontWeight','bold','Color',grey_connections);
text(212.5, 1270,'C3','FontSize',fontsize+6,'FontName','Avenir','FontAngle','italic',...
    'FontWeight','bold','Color',grey_connections);
text(207.3, 1340,'C2','FontSize',fontsize+6,'FontName','Avenir','FontAngle','italic',...
    'FontWeight','bold','Color',grey_connections);

yyaxis right
hold on;
xlim([155 216]); ylim([0.0 1.1]);
set(gca,'ycolor',[0    0.4470    0.7410]); set(gca,'ycolor','k','FontSize',fontsize);
set(gca, 'SortMethod', 'depth','ytick',[])
lgd=legend('Formation','{\it v_{crit}}','HF event','Location','NorthEast','Orientation','vertical');
set(lgd,'Box','on','FontName','Avenir','FontSize',fontsize-1,'LineWidth',0.5,...
    'Position',[0.695 0.915 0.01 0.01]);
grid on; box on
ax=gca; ax.Layer = 'top';

%%%%%%%%%% log cluster sizes + how many days it takes to make a cluster that large %%%%%%%%%%%%%%%%%%%%%%
[~, ~, contributing_rows, cluster_num_lakes,...
    cluster_num_timepoints] = count_simultaneous_drainage(days_HF_drainages2);
cluster_time_window_S2_2022 = cluster_num_timepoints; % [ days ]
cluster_size_S2_2022 = cluster_num_lakes; % [ number of lakes ]
% keep unique clusters: must have same exact contributing_rows
strCellArray = cellfun(@(x) mat2str(sort(x)), contributing_rows, 'UniformOutput', false);
% Find unique sorted rows
[uniqueStr, ~, ~] = unique(strCellArray);
% Initialize a vector to store row numbers of unique rows
uniqueRowNumbers = [];
% Loop through unique rows to save their original row numbers
for i = 1:numel(uniqueStr)
    rowNumbers = find(strcmp(strCellArray, uniqueStr{i}));
    uniqueRowNumbers = [uniqueRowNumbers rowNumbers(1)];
end
% Sort uniqueRowNumbers to maintain the order of appearance
uniqueRowNumbers = sort(uniqueRowNumbers);
% Get unique rows in cluster time and sizes
unique_HF_cluster_S2_2022(:,1) = cluster_time_window_S2_2022(uniqueRowNumbers);
unique_HF_cluster_S2_2022(:,2) = cluster_size_S2_2022(uniqueRowNumbers);

for i=1:1:length(m)
    for j=1:1:length(t)
    PDF(i,j) = (((lambda.*t(j)).^(m(i))).*(exp(-1.*lambda.*t(j))))./...
        (factorial(m(i))); % make PDF for all m and t 
    end
end
% probability of observing one or more events less than or equal to size m = 
% 1 minus the probability of observing zero events + an event of size 1 +
% + an event of size 2 + ... + an event of size (m-1) + an event of size m:
CDF_2 = 1-cumsum(PDF,1); % CDF
% find values at a given \alpha threshold: (it's about ~1.23 * t * \lambda)
for i=1:1:length(t)
    [sig_HF_cluster(i,1)] = find(CDF_2(:,i)<=alpha,1)-1;
end
% clusters observed in 2022, as a f(cluster size, cluster probability),
% given the time window observed
clear cluster_probability_S2_2022
for i=1:1:length(unique_HF_cluster_S2_2022(:,1))
    cluster_probability_S2_2022(i,1) = CDF_2(unique_HF_cluster_S2_2022(i,2)+1,... % num of events
                                             unique_HF_cluster_S2_2022(i,1)); % time window
end

%%%%%%%%%% little cluster plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(axe2); 
hold on;
plot(0:1:15,(0:1:15).*lambda,'-','Color',[0.8500    0.3250    0.0980],'LineWidth',1.2);
% significant region
a1 = area(0:1:length(sig_HF_cluster),vertcat(0,sig_HF_cluster));
a1.FaceColor = lake_blue(1:3); a1.FaceAlpha = 0.15;
a1.EdgeColor = [0.8500    0.3250    0.0980]; a1.BaseLine.LineStyle = 'none';
a1.BaseValue = 15; a1.LineWidth = 1.9; a1.LineStyle = ':';
% all clusters
plot(unique_HF_cluster_S2_2022(:,1),unique_HF_cluster_S2_2022(:,2),'dk',... 
    'MarkerFaceColor','none','MarkerSize',5); 
% significant clusters
for i=1:length(unique_HF_cluster_S2_2022(:,1))
    if cluster_probability_S2_2022(i) <= alpha
        plot(unique_HF_cluster_S2_2022(i,1),unique_HF_cluster_S2_2022(i,2),'kp',...
        'MarkerFaceColor',[247 252 185]./255,'MarkerSize',10);
    else
    end
end
% C2 northern 
plot(4,3,'dk','MarkerFaceColor','none','MarkerSize',5); % C2.North 
quiver(4, 5.8, 0, -2.8,'LineWidth',1.5,'Color',lake_blue(1:3),'MaxHeadSize',0.5);
text(4,2.5,'C2.N','FontSize',fontsize+1,'FontName','Avenir','FontAngle','italic',...
    'FontWeight','bold','Color',grey_connections);
% C2 southern
%plot(2,3,'dk','MarkerFaceColor','none','MarkerSize',4); % C2.South
quiver(4, 5.8, -2.1, -2.9,'LineWidth',1.5,'Color',lake_blue(1:3),'MaxHeadSize',0.5);
text(1.5,2.4,'C2.S','FontSize',fontsize+1,'FontName','Avenir','FontAngle','italic',...
    'FontWeight','bold','Color',grey_connections);

text(0.3,6.5,'b','FontName','Arial','FontWeight','bold','FontSize',14)
text(0.25,5.55,'C1','FontSize',fontsize+4,'FontName','Avenir','FontAngle','italic',...
    'FontWeight','bold','Color',grey_connections);
text(1.4,4.6,'C3','FontSize',fontsize+4,'FontName','Avenir','FontAngle','italic',...
    'FontWeight','bold','Color',grey_connections);
text(3,6.5,'C2','FontSize',fontsize+4,'FontName','Avenir','FontAngle','italic',...
    'FontWeight','bold','Color',grey_connections);
text(4.3,3.55,sprintf('\\lambda = %3.2f',lambda),'FontName','Arial','FontSize',fontsize,'rotation',34,...
    'Color',[0.8500    0.3250    0.0980],'FontWeight','bold');
ylabel('Number of HF events, {\itm}'); xlabel('Time window, {\itw} [ days ]'); 
ylim([0 7]); xlim([0 6]); grid on; box on;
set(gca,'LineWidth',1.1,'tickdir','in','TickLength',[0 0],'FontSize',fontsize-1,'YAxisLocation','right','ytick',0:1:6,'xtick',0:1:10);


%% %%%%%%%% drainage mechanism by elevation plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(axe3)
yyaxis left
set(gca,'YtickLabel',[],'YColor','k')
set(gca,'ytick',700:100:1800,'xtick',0:25:100,'xTickLabel',0:25:100)

yyaxis right
hold on;
test = struct2table(counting_elev_lakes_ROI_22);
b = barh(test.elev_lower_bound+50, test.percent_form, 'stacked','EdgeColor',[0.6 0.6 0.6],'BarWidth',1);
b(1).FaceColor = goldenrod;
b(2).FaceColor = moulin_green;
b(3).FaceColor = lake_blue(1:3);
b(4).FaceColor = plum(1:3);
set(gca,'YColor','k')
hold all

text(50, 1740, [{'No lakes'};{'formed here.'}],'FontName','Avenir','FontAngle','italic',...
    'FontSize',fontsize,'HorizontalAlignment','center');

for i=1:length(test.elev_lower_bound)-2
    number_plot = double(test.number_form(i));
    elev_plot = double(test.elev_lower_bound(i))+50;
    text(110, elev_plot, num2str(number_plot),'FontSize',fontsize-1,'FontName','AvenirBold',...
        'HorizontalAlignment','center')
end

lgd = legend('Hydro-fracture','Moulin','Overspill','No-exit, frozen');
set(lgd,'Position',[0.585 0.965 0.01 0.01],'NumColumns',2,'Box','on',...
    'FontName','Avenir','FontSize',fontsize-1,'LineWidth',0.5)

text(5, 1650, 'c','FontName','Arial','FontWeight','bold','FontSize',fontsize+4,'Color','w');
xlabel('Lake Fate [ % ]')
xlim([0 100]); ylim([700 1800]); 
set(gca,'ytick',700:100:1800,'xtick',0:25:100,'xTickLabel',0:25:100,'yTickLabel',[])
set(gca,'YColor','k')
box on; set(gca,'LineWidth',1.1,'FontName','Avenir','FontSize',fontsize-2);
ax=gca; ax.Layer = 'top';


%% %%%%%%%%%%%%%%%%% 2023 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all 

% CLASSIFIER
% 1 = HF; 2 = moulin; 3 = overspill; 4 = crevasses; 5 = frozen. 9 = throw out/inconclusive.
load('daily_lake_HFpossible_classifier_S1S2MO_2023_251003.mat'); % lakes large enough to HF
daily_lake_classifier = daily_lake_HFpossible_classifier;
 
% lakes large enough to HF and in nevis domain:
% load lake centre points, boundaries from FASTER for 2023
%  boundaries and correct dates:
load('environs_lakes_2023B_S1S2MO_251003.mat') %  boundaries and correct dates
lake_typeing = environs_lakes.laketypeing_dates; % has NaN  outside the nevis area
lake_typeing(environs_lakes.lake_HF_possible_days==0,7) = NaN; % 0 days of HF possible
lake_elev = environs_lakes.S; % lake elevation from BedMachine

% load available image dates for FASTER for 2023
load('all_JD_nums_S1S2_2023.mat'); %all julian day numbers
load('imageID_2023_250523.mat'); % JD numbers of S2 images available for each specific lake

% ROI that's relevant to the GNSS array after FASTER locations update (2025-April):
ROI_relevant_xv = [-15 60 104 104 15 -10 -33 -33 -15]; % km in nevis-model space
ROI_relevant_yv = [15 4 -8 -50 -50 -41 -33 -1 15]; % km in nevis-model space
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

% lake_typeing inside the ROI (for elevation-bin analysis 25-02-12)
lake_typeing_all = lake_typeing;
% count up lakes by elevation bands
elev_band_interval = 100; 
elev_bands = 700:elev_band_interval:1900; % 100 m a.s.l.
for i=1:1:length(elev_bands)-1
    % find lakes within this elevation band
    counting_elev_lakes(i).id = find(lake_elev>elev_bands(i) & lake_elev<=elev_bands(i+1));
    % count how many HF/moulin/overspill/frozen events there are
    counting_elev_lakes(i).mechanism(1,1) = length(find(lake_typeing_all(counting_elev_lakes(i).id,7)==1));
    counting_elev_lakes(i).mechanism(1,2) = length(find(lake_typeing_all(counting_elev_lakes(i).id,7)==2));
    counting_elev_lakes(i).mechanism(1,3) = length(find(lake_typeing_all(counting_elev_lakes(i).id,7)==3));
    counting_elev_lakes(i).mechanism(1,4) = length(find(lake_typeing_all(counting_elev_lakes(i).id,7)==5));
    % count how many lakes there are
    counting_elev_lakes(i).number_form = sum(counting_elev_lakes(i).mechanism(1,1:4));
    % percentage of HF out of all lakes (aka lakes that form)
    counting_elev_lakes(i).percent_form = 100.*(counting_elev_lakes(i).mechanism./counting_elev_lakes(i).number_form);
    % count how many lakes drain
    counting_elev_lakes(i).number_drain = sum(counting_elev_lakes(i).mechanism(1,1:3));
    % percentage of HF out of drainage events (aka lakes that drain)
    counting_elev_lakes(i).percent_drain = 100.*(counting_elev_lakes(i).mechanism(1:3)./counting_elev_lakes(i).number_drain);
    % save elevation band lower bound
    counting_elev_lakes(i).elev_lower_bound = elev_bands(i);
end
counting_elev_lakes_ROI_23 = counting_elev_lakes;

% racmo at time of HF drainage
load cumulative_runoff_last_day_to_HF_2023.mat

% preliminaries for homo poisson
alpha = 0.05; % significance threshold
t = 1:1:60; % time window [ days ]
m = 0:1:250; % number of events [ cluster size ]

% HF 2023 meltseason 
id_HF_drainages = find(lake_typeing(:,7)==1); % HF drainages row IDs [ HF == 1 ]
save id_HF_drainages_2023_250415.mat id_HF_drainages
drainages = length(id_HF_drainages)-1 % number of HF drainages in 2023 (take out one at 226–229)
days_HF_drainages = daily_lake_classifier(id_HF_drainages,:); % HF drainages timeseries
days_HF_drainages(days_HF_drainages~=2) = NaN; % NaN if not in drainage window
days_HF_drainages(days_HF_drainages==2) = 1; % drainage window day ---> 1

% % get the in-between days, and not the initial day, --> 1
days_HF_drainages2 = NaN(length(id_HF_drainages),365); % prefill with NaNs
for i=1:1:length(id_HF_drainages)
    % if just one day draining window: it's that day
    if environs_lakes.laketypeing_dates(id_HF_drainages(i),3) == environs_lakes.laketypeing_dates(id_HF_drainages(i),4)
        days_HF_drainages2(i,environs_lakes.laketypeing_dates(id_HF_drainages(i),4)) = 1;
    elseif environs_lakes.laketypeing_dates(id_HF_drainages(i),3)+1 == environs_lakes.laketypeing_dates(id_HF_drainages(i),4)
        days_HF_drainages2(i,environs_lakes.laketypeing_dates(id_HF_drainages(i),4)) = 1;
    else 
        days_HF_drainages2(i,(environs_lakes.laketypeing_dates(id_HF_drainages(i),3)+1):environs_lakes.laketypeing_dates(id_HF_drainages(i),4)) = 1;
   end
end
days_HF_drainages = days_HF_drainages2;

% get the in-between days, and not the initial day, --> 1
%days_HF_drainages = diff(days_HF_drainages,1,2) + 1;  % diffs along rows
%figure(1); imagesc(days_HF_drainages); colorbar; xlim([160 280]); shg
% drainage rates
drainages_start = 184; % [ DOY ]
drainages_stop = 206; % [ DOY ]
days = drainages_stop-drainages_start; % day window HF drainages occur on in 2023
lambda = (drainages)./days; 
% calculate number of HF lakes draining on all days
days_HF_drainages_sum = nansum(days_HF_drainages,1);
% calculate error in each lake-drainage date (this is the window size)
days_HF_drainages_window = nansum(days_HF_drainages,2);
average_HF_window = nanmean(days_HF_drainages_window);
% lake elevations
elevations_HF = lake_elev(id_HF_drainages);

%% %%%%%%%%% v_crit and lake formation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get v_crit and lake-formation time windows
days_HF_form = daily_lake_classifier(id_HF_drainages,:); % HF drainages timeseries
for i=1:length(days_HF_form(:,1))
    % whats the first image we have of this lake?
    day_lake_form_last(i,1) = lake_typeing(id_HF_drainages(i,1),2); % image date of lake first formed
    index1 = find(all_JD_nums(:,1)==day_lake_form_last(i,1)); % Sentinel image index on that day
    if index1>1
    index2 = find(imageID_2023.lake(id_HF_drainages(i),1).plot_IDs>=index1,1,'first'); % Sentinel image on day lake first formed (lake specific)
    % what's the last image available before we see a lake?  (date of last dry basin) 
    index3 = imageID_2023.lake(id_HF_drainages(i),1).plot_IDs(index2-1); % image index before lake first formed (lake specific)
    day_lake_form_first(i,1) = all_JD_nums(index3,1)+1; % date one day after image of dry lake basin (lake specific)
    error_day_lake_form(i,1) = day_lake_form_last(i,1)-day_lake_form_first(i,1); % days of spread (count 24-hr periods)
    else
    day_lake_form_first(i,1) = NaN; % don't plot a formation range if there is no prior image (this is just a problem for 1 HF lake)
    end
    % v_crit model estimate
    day_HF_possible(i,1) = find(days_HF_form(i,:)==1,1,'first'); % first day of lake being HFpossible
end
% for plotting
formation_average = (day_lake_form_first+day_lake_form_last)/2;
formation_error = day_lake_form_last-formation_average;

%% plotting  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axe4 = axes('Position',[0.085 0.042 0.705 0.45],'Box','on','NextPlot','add');
axe5 = axes('Position',[0.095 0.317 0.23 0.17],'Box','on','NextPlot','add');
axe6 = axes('Position',[0.80 0.042 0.17 0.45],'Box','on','NextPlot','add');

fontsize = 9; 
% more colors
ice_sheet_contours = [0.7 0.7 0.7];
plum = [211, 160, 211, 100]./255;
lake_blue = [0    0.4470    0.7410   0.08]; % fourth value is FaceAlpha built in!
goldenrod = [1.0 0.8 0.0];
icy_plum = [203, 183, 227]./255;
moulin_green = [0.4 0.8 0.8];
lake_colors = vertcat(goldenrod, moulin_green, lake_blue(1:3), icy_plum);

axes(axe4); hold on
yyaxis left

% RACMO : passing melt thresholds at elevation bands
load('runoff_elev_2023_300m.mat')
melt_threshold = [50, 100, 300, 500, 700, 900, 1100, 1300]; % 200 mm = 0.2 m
melt_threshold = [50, 250, 450, 650, 850, 1050, 1250]; 
for i=1:1:length(runoff_elev_2023_300m.elev_bands)-1
    for j=1:length(melt_threshold)
    if runoff_elev_2023_300m.cumulative(i,end) >= melt_threshold(j)
    DOY_threshold(i,j) = find(runoff_elev_2023_300m.cumulative(i,:)>melt_threshold(j),1,'first'); % first day above 0.25 m
    else
    DOY_threshold(i,j) = 1e9;
    end
   end
end
melt_colors = [[0.9 0.9 0.9],;
       [0.85 0.85 0.85],;
       [0.8 0.8 0.8],;
       [0.75 0.75 0.75],;
       [0.7 0.7 0.7],;
       [0.65 0.65 0.65],;
       [0.5 0.5 0.5]];
area(DOY_threshold(:,1),runoff_elev_2023_300m.elev_bands(1:end-1),'EdgeColor','none','FaceColor',melt_colors(1,:)); 
area(DOY_threshold(:,2),runoff_elev_2023_300m.elev_bands(1:end-1),'EdgeColor','none','FaceColor',melt_colors(2,:)); 
area(DOY_threshold(:,3),runoff_elev_2023_300m.elev_bands(1:end-1),'EdgeColor','none','FaceColor',melt_colors(3,:)); 
area(DOY_threshold(:,4),runoff_elev_2023_300m.elev_bands(1:end-1),'EdgeColor','none','FaceColor',melt_colors(4,:)); 
area(DOY_threshold(:,5),runoff_elev_2023_300m.elev_bands(1:end-1),'EdgeColor','none','FaceColor',melt_colors(5,:)); 
area(DOY_threshold(:,6),runoff_elev_2023_300m.elev_bands(1:end-1),'EdgeColor','none','FaceColor',melt_colors(6,:));
area(DOY_threshold(:,7),runoff_elev_2023_300m.elev_bands(1:end-1),'EdgeColor','none','FaceColor',melt_colors(7,:));
% labels vertical
text(214.75,1670,{'Runoff:'},'FontSize',fontsize+1,'FontName','Avenir','Rotation',0,'HorizontalAlignment','right');
text(214.15,1540,sprintf('%d cm',melt_threshold(1)./10),'FontName','Avenir','Rotation',90);
text(214.15,1400,sprintf('%d cm',melt_threshold(2)./10),'FontName','Avenir','Rotation',90);
text(214.15,1270,sprintf('%d cm',melt_threshold(3)./10),'FontName','Avenir','Rotation',90);
text(214.15,1140,sprintf('%d cm',melt_threshold(4)./10),'FontName','Avenir','Rotation',90);
text(214.15,965,sprintf('%d cm',melt_threshold(5)./10),'FontName','Avenir','Rotation',90);
text(214.15,765,sprintf('%d cm',melt_threshold(6)./10),'FontName','Avenir','Rotation',90);
    
% vcrit, formation dates
errorbar(formation_average, elevations_HF,[],[],formation_error,formation_error,'k.',...
    'LineWidth',0.7,'Marker','none','CapSize',4); % lake formation
plot(day_HF_possible,elevations_HF,'dk',...
    'markerFaceColor','w','MarkerSize',5); % v_crit

% HF clusters
rectangle('Position',[190.55 820 2.8 340],'EdgeColor',lake_blue(1:3),'FaceColor','none','LineWidth',1.1);
rectangle('Position',[200.5 1300 1 240],'EdgeColor',lake_blue(1:3),'FaceColor','none','LineWidth',1.1);
plot([194.5 194.5 200.5 200.5 201.5 201.5 194.5],...
    [1160 1325 1325 1300 1300 1160 1160],'-','Color',lake_blue(1:3),'LineWidth',1.1)
for i=1:1:length(days_HF_drainages(:,1))
    plot(1:1:365,days_HF_drainages(i,:)+elevations_HF(i)-days_HF_drainages(i,:),'-k',...
        'LineWidth',0.8); 
end
for i=1:1:length(days_HF_drainages(:,1))
    plot(1:1:365,days_HF_drainages(i,:)+elevations_HF(i)-days_HF_drainages(i,:),'dk',...
        'markerFaceColor',goldenrod,'MarkerSize',5); 
end
xlim([155 216]); ylim([700 1800]);
set(gca,'ytick',500:100:1800);
set(gca,'LineWidth',1.1,'tickdir','in','TickLength',[0 0],'FontSize',fontsize,'FontName','Avenir');
set(gca,'ycolor','k');
ylabel('Ice-sheet Surface Elevation [ m a.s.l. ]','FontSize',fontsize);
xlabel('Day of Year, 2023','FontSize',fontsize)
%title('Hydro-fracture drainages constrained by 2023 Sentinel-2 images','FontSize',11)
set(gca, 'SortMethod', 'depth')
text(148.5, 1800, 'd','FontName','Arial','FontWeight','bold','FontSize',fontsize+4);
% % instrumented lakes
% % 86, 81
% text(find(days_HF_drainages(9,:)==1,1,'first')+3.5, elevations_HF(9)-7, 'L1A','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10);
% text(find(days_HF_drainages(8,:)==1,1,'first')+3.5, elevations_HF(8)+17, 'L1B','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10);
% % % 142, 132, 150
% text(find(days_HF_drainages(13,:)==1,1,'first')+5.5, elevations_HF(13), 'L2A','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10);
% text(find(days_HF_drainages(12,:)==1,1,'first')+5.5, elevations_HF(12)-5, 'L2B','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10);
% text(find(days_HF_drainages(15,:)==1,1,'first')-1.75, elevations_HF(15), 'L150','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10);
% % 208, 201, 228, 238, 273, 275
% text(201.75, elevations_HF(25)+8, 'L3A','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10);
% text(201.75, elevations_HF(24)-11, 'L3B','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10);
% text(201.75, elevations_HF(27)-18, 'L3C','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10);
% text(201.75, elevations_HF(27)+8, 'L3D','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10);
% text(201.75, elevations_HF(28)-2, 'L3E','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10);
% text(201.75, elevations_HF(29)+5, 'L3F','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10);
 
% % other lakes
% text(find(days_HF_drainages(1,:)==1,1,'first')+0.5, elevations_HF(1), 'L8','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% % text(find(days_HF_drainages(2,:)==1,1,'first')+5.25+0.5, elevations_HF(2), 'L12','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% % text(find(days_HF_drainages(3,:)==1,1,'first')+5.25+0.5, elevations_HF(3), 'L28','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages(2,:)==1,1,'first')-0.2+0.5, elevations_HF(2), 'L38','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages(3,:)==1,1,'first')+5, elevations_HF(3)-10, 'L44','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% % text(find(days_HF_drainages(6,:)==1,1,'first')+0.5, elevations_HF(6), 'L47','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages(4,:)==1,1,'first')+0.75, elevations_HF(4)+14, 'L56','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages(5,:)==1,1,'first')-0.2+0.5, elevations_HF(5), 'L59','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages(6,:)==1,1,'first')+5, elevations_HF(6)+7, 'L66','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages(7,:)==1,1,'first')+4.5+0.5, elevations_HF(7), 'L75','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% % text(find(days_HF_drainages(12,:)==1,1,'first')+0.5, elevations_HF(12), 'L85','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% % text(find(days_HF_drainages(14,:)==1,1,'first')+0.5, elevations_HF(14), 'L88','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages(10,:)==1,1,'first')-0.25+0.5, elevations_HF(10)+5, 'L103','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% % text(find(days_HF_drainages(16,:)==1,1,'first')-0.25+0.5, elevations_HF(16)-20, 'L119','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% % text(find(days_HF_drainages(17,:)==1,1,'first')-0.25+0.5, elevations_HF(17)-8, 'L123','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages(11,:)==1,1,'first')+0.5, elevations_HF(11), 'L131','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% %text(find(days_HF_drainages(21,:)==1,1,'first')+0.5, elevations_HF(21)-15, 'L146','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages(16,:)==1,1,'first')+9+0.5, elevations_HF(16)+5, 'L152','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages(17,:)==1,1,'first')+0.25, elevations_HF(17)+15, 'L154','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages(18,:)==1,1,'first')+7+0.5, elevations_HF(18)+35, 'L162','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% % text(find(days_HF_drainages(26,:)==1,1,'first')+4+0.5, elevations_HF(26)-4, 'L166','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages(19,:)==1,1,'first')+0.25, elevations_HF(19)-10, 'L171','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages(20,:)==1,1,'first')+0.5, elevations_HF(20)-40, 'L179','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages(21,:)==1,1,'first')+2+0.5, elevations_HF(21)+25, 'L181','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% % text(find(days_HF_drainages(30,:)==1,1,'first')+2+0.5, elevations_HF(30)+45, 'L187','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% % text(find(days_HF_drainages(31,:)==1,1,'first')+5+0.5, elevations_HF(31), 'L189','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% % text(find(days_HF_drainages(32,:)==1,1,'first')+4+0.5, elevations_HF(32)+10, 'L191','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages(22,:)==1,1,'first')+7+0.5, elevations_HF(22), 'L194','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% text(find(days_HF_drainages(23,:)==1,1,'first')+0.5, elevations_HF(23)-20, 'L196','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% % text(find(days_HF_drainages(38,:)==1,1,'first')+0.5, elevations_HF(38)-1, 'L230','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% % text(find(days_HF_drainages(39,:)==1,1,'first')+0.5, elevations_HF(39)+24, 'L235','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% % text(find(days_HF_drainages(41,:)==1,1,'first')+7+0.5, elevations_HF(41), 'L268','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% % text(find(days_HF_drainages(44,:)==1,1,'first')+5+0.5, elevations_HF(44)-8, 'L288','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% % text(find(days_HF_drainages(45,:)==1,1,'first')+7+0.5, elevations_HF(45), 'L299','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% % text(find(days_HF_drainages(46,:)==1,1,'first')+5+0.5, elevations_HF(46)+15, 'L310','FontName','Avenir','FontWeight','bold','FontAngle','italic','FontSize',10,'HorizontalAlignment','right');
% cluster labels
grey_connections = lake_blue(1:3).*1.3; %[0.7 0.7 0.7];
text(190.65, 785,'C4','FontSize',fontsize+6,'FontName','Avenir','FontAngle','italic',...
    'FontWeight','bold','Color',grey_connections);
text(199.9, 1575,'C6','FontSize',fontsize+6,'FontName','Avenir','FontAngle','italic',...
    'FontWeight','bold','Color',grey_connections);
text(196.3, 1360,'C5','FontSize',fontsize+6,'FontName','Avenir','FontAngle','italic',...
    'FontWeight','bold','Color',grey_connections);

yyaxis right
hold on;
set(gca,'ycolor',[0    0.4470    0.7410]); set(gca,'ycolor','k','FontSize',fontsize);
set(gca, 'SortMethod', 'depth','ytick',[])
grid on; box on; 
ax=gca; ax.Layer = 'top';

[num_drainage_timepoints, drainage_indices, contributing_rows, cluster_num_lakes,...
    cluster_num_timepoints] = count_simultaneous_drainage(days_HF_drainages);
cluster_time_window_S2_2022 = cluster_num_timepoints; % [ days ]
cluster_size_S2_2022 = cluster_num_lakes; % [ number of lakes ]
% keep unique clusters: must have same exact contributing_rows
strCellArray = cellfun(@(x) mat2str(sort(x)), contributing_rows, 'UniformOutput', false);
% Find unique sorted rows
[uniqueStr, ~, idx] = unique(strCellArray);
% Initialize a vector to store row numbers of unique rows
uniqueRowNumbers = [];
% Loop through unique rows to save their original row numbers
for i = 1:numel(uniqueStr)
    rowNumbers = find(strcmp(strCellArray, uniqueStr{i}));
    uniqueRowNumbers = [uniqueRowNumbers rowNumbers(1)];
end
% Sort uniqueRowNumbers to maintain the order of appearance
uniqueRowNumbers = sort(uniqueRowNumbers);
% Get unique rows in cluster time and sizes
unique_HF_cluster_S2_2022(:,1) = cluster_time_window_S2_2022(uniqueRowNumbers);
unique_HF_cluster_S2_2022(:,2) = cluster_size_S2_2022(uniqueRowNumbers);

for i=1:1:length(m)
    for j=1:1:length(t)
    PDF(i,j) = (((lambda.*t(j)).^(m(i))).*(exp(-1.*lambda.*t(j))))./...
        (factorial(m(i)));
    end
end
% probability of observing one or more events less than or equal to size m = 
% 1 minus the probability of observing zero events + an event of size 1 +
% + an event of size 2 + ... + an event of size (m-1) + an event of size m:
CDF_2 = 1-cumsum(PDF,1);
% find values at a given \alpha threshold: (~1.23 * t * \lambda)
for i=1:1:length(t)
    [sig_HF_cluster(i,1)] = find(CDF_2(:,i)<=alpha,1)-1;
end
% clusters observed in 2022, as a f(cluster size, cluster probability),
% given the time window observed
clear cluster_probability_S2_2022
for i=1:1:length(unique_HF_cluster_S2_2022(:,1))
    cluster_probability_S2_2022(i,1) = CDF_2(unique_HF_cluster_S2_2022(i,2)+1,... % num of events
                                             unique_HF_cluster_S2_2022(i,1)); % time window
end

%%%%%%%%%%%%%%%%%%%% little clusters plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(axe5); 
hold on;
plot(0:1:15,(0:1:15).*lambda,'-','Color',[0.8500    0.3250    0.0980],'LineWidth',1.2);
% significant area
a1 = area(0:1:length(sig_HF_cluster),vertcat(0,sig_HF_cluster));
a1.FaceColor = [0.8500    0.3250    0.0980]; a1.FaceAlpha = 0.25;
a1.FaceColor = lake_blue(1:3); a1.FaceAlpha = 0.15;
a1.EdgeColor = [0.8500    0.3250    0.0980]; a1.BaseLine.LineStyle = 'none';
a1.BaseValue = 50; a1.LineWidth = 1.9; a1.LineStyle = ':';

plot(unique_HF_cluster_S2_2022(:,1),unique_HF_cluster_S2_2022(:,2),'dk',... 
    'MarkerFaceColor','none','MarkerSize',5); 
% significant clusters
for i=1:length(unique_HF_cluster_S2_2022(:,1))
    if cluster_probability_S2_2022(i) <= alpha
        plot(unique_HF_cluster_S2_2022(i,1),unique_HF_cluster_S2_2022(i,2),'kp',...
        'MarkerFaceColor',[247 252 185]./255,'MarkerSize',10); %
    else
    end
end
% C5 northern 
%plot(5,10,'dk','MarkerFaceColor','none','MarkerSize',4); % C5.North and C5.South
quiver(7, 15.4, -2.1, -5.4,'LineWidth',1.5,'Color',lake_blue(1:3),'MaxHeadSize',0.4);
text(4.25,8.75,'C5.N','FontSize',fontsize+1,'FontName','Avenir','FontAngle','italic',...
    'FontWeight','bold','Color',grey_connections);
% C5 southern
plot(7,11,'dk','MarkerFaceColor','none','MarkerSize',5); % C5.North and C5.South
quiver(7, 15.4, 0, -4.2,'LineWidth',1.5,'Color',lake_blue(1:3),'MaxHeadSize',0.5);
text(6,9.9,'C5.S','FontSize',fontsize+1,'FontName','Avenir','FontAngle','italic',...
    'FontWeight','bold','Color',grey_connections);

text(0.3,15.25,'e','FontName','Arial','FontWeight','bold','FontSize',fontsize+4)
text(2.25,9.5,'C4','FontSize',fontsize+4,'FontName','Avenir','FontAngle','italic',...
    'FontWeight','bold','Color',grey_connections);
text(5.5,15.8,'C5','FontSize',fontsize+4,'FontName','Avenir','FontAngle','italic',...
    'FontWeight','bold','Color',grey_connections);
text(0.5,7.5,'C6','FontSize',fontsize+4,'FontName','Avenir','FontAngle','italic',...
    'FontWeight','bold','Color',grey_connections);
text(2.6,4.6,sprintf('\\lambda = %3.2f',lambda),'FontName','Arial','FontSize',fontsize,'rotation',36,...
    'Color',[0.8500    0.3250    0.0980],'FontWeight','bold');
ylabel('Number of HF events, {\itm}'); xlabel('Time window, {\itw} [ days ]'); 
ylim([0 17]); xlim([0 8]); grid on; box on;
set(gca,'LineWidth',1.1,'tickdir','in','TickLength',[0 0],'FontSize',fontsize-1,...
    'ytick',0:2:16,'xtick',0:1:10,'YAxisLocation','right');


%%%%%%%%%% drainage mechanism by elevation plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(axe6)
yyaxis left
set(gca,'YtickLabel',[],'YColor','k')
set(gca,'ytick',700:100:1800,'xtick',0:25:100,'xTickLabel',0:25:100)

yyaxis right
hold on;
test = struct2table(counting_elev_lakes_ROI_23);
b = barh(test.elev_lower_bound+50, test.percent_form, 'stacked','EdgeColor',[0.6 0.6 0.6],'BarWidth',1);
b(1).FaceColor = goldenrod;
b(2).FaceColor = moulin_green;
b(3).FaceColor = lake_blue(1:3);
b(4).FaceColor = plum(1:3);
set(gca,'YColor','k')
hold all

% number of lakes 
for i=1:length(test.elev_lower_bound)-1
    number_plot = double(test.number_form(i));
    elev_plot = double(test.elev_lower_bound(i))+50;
    text(110, elev_plot, num2str(number_plot),'FontSize',fontsize-1,'FontName','AvenirBold',...
        'HorizontalAlignment','center')
end

text(5, 1750, 'f','FontName','Arial','FontWeight','bold','FontSize',fontsize+4,'Color','w');
xlabel('Lake Fate [ % ]','FontSize',fontsize)
xlim([0 100]); ylim([700 1800]); 
set(gca,'ytick',700:100:1800,'xtick',0:25:100,'xTickLabel',0:25:100,'yTickLabel',[])
set(gca,'YColor','k')
box on; set(gca,'LineWidth',1.1,'FontName','Avenir','FontSize',fontsize-2);
ax=gca; ax.Layer = 'top';

%% print figure
print(gcf,'-dpng','-r600',sprintf('../paperfigs/paperfig3_homopo_vcrit_HF_260128_repository.png')); 
