%% clusters - homogeneous Poisson process for OVERSPILL (OS) 
% 24-08-23 LAS -- slim panels down for paper figures :D
% 24-09-03 LAS -- time of drainage plotted as DOY values of 24-hr periods
% 24-09-03 LAS -- Change question to ask what's stat sig in the drainage
% catchment? (i.e., just the region of the ROI that's relevant to the GNSS
% array)
% 24-10-28 LAS -- racmo accumulated over time
% 25-04-16 LAS -- update with FASTER lake boundaries and locations
% 25-05-28 LAS -- v_crit and lake formation (in DOY)
% 26-01-28 LAS -- repository checks

clear all; close all;

%% CLASSIFIER 2022
% 1 = HF; 2 = moulin; 3 = overspill; 4 = crevasses; 5 = frozen. 9 = throw out/inconclusive.
load('daily_lake_HFpossible_classifier_S1S2_2022_240829.mat'); % S1+S2 images

% lakes large enough to HF and in nevis domain:
% load lake boundaries, centre points, drainage dates from FASTER for 2022
load('environs_lakes_2022B_250416.mat') % drainage dates & boundaries
environs_lakes_2022 = environs_lakes; 
lake_typeing = environs_lakes.laketypeing_dates; % has NaN outside the nevis area
lake_typeing(environs_lakes.lake_HF_possible_days==0,7) = NaN; % 0 days of HF possible
lake_elev = environs_lakes.S; % lake elevation from BedMachine

% load available image dates for FASTER for 2022
load('all_JD_nums_2022.mat'); % all julian day numbers of S2 images
load('imageID_2022_250523.mat'); % JD numbers of S2 images available for each specific lake

% ROI that's relevant to the GNSS array, after FASTER locations update (2025-April):
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
% NaN from the classifier
environs_lakes.laketypeing_dates(180,1:7) = [NaN, NaN, NaN, NaN, NaN, NaN, NaN];
environs_lakes.laketypeing_dates(208,1:7) = [NaN, NaN, NaN, NaN, NaN, NaN, NaN];
% redefine lake_typeing to just include lakes within nevis domain
lake_typeing = environs_lakes.laketypeing_dates; % after throwing out lakes on the nevis domain border

% polygon ROI check
% figure; plot(ROI_relevant_xv,ROI_relevant_yv)
% hold on; plot(environs_lakes.X_km,environs_lakes.Y_km,'o')
% %pause; close all

% racmo at time of HF
load cumulative_runoff_last_day_to_HF_2022.mat

% preliminaries for homo poisson
alpha = 0.05; % significance threshold
t = 1:1:70; % time window [ days ]
m = 0:1:250; % number of events [ cluster size ]
% drainage dates for OVERSPILL 2022
id_OS_drainages = find(lake_typeing(:,7)==3); % OVERSPILL drainages row IDs [ OS == 3 ]
days_OS_drainages = daily_lake_classifier(id_OS_drainages,:); % OS drainages timeseries
days_OS_drainages(days_OS_drainages~=2) = NaN; % NaN if not in filling or drainage window
% change 2's to 1's
days_OS_drainages(days_OS_drainages==2) = 1; % drainage window day ---> 1's
% remove rows that are all NaNs
days_OS_drainages = days_OS_drainages(~all(isnan(days_OS_drainages), 2), :);
%figure(5); imagesc(days_OS_drainages); colorbar; xlim([160 260]); title('OS start stop'); shg

% % get the in-between days, and not the initial day, --> 1
days_OS_drainages2 = NaN(length(id_OS_drainages),365); % prefill with NaNs
for i=1:1:length(id_OS_drainages)
    % if just one day draining window: it's that day
    if environs_lakes.laketypeing_dates(id_OS_drainages(i),3) == environs_lakes.laketypeing_dates(id_OS_drainages(i),4)
        days_OS_drainages2(i,environs_lakes.laketypeing_dates(id_OS_drainages(i),4)) = 1;
    elseif environs_lakes.laketypeing_dates(id_OS_drainages(i),3)+1 == environs_lakes.laketypeing_dates(id_OS_drainages(i),4)
        days_OS_drainages2(i,environs_lakes.laketypeing_dates(id_OS_drainages(i),4)) = 1;
    else 
        days_OS_drainages2(i,(environs_lakes.laketypeing_dates(id_OS_drainages(i),3)+1):environs_lakes.laketypeing_dates(id_OS_drainages(i),4)) = 1;
   end
end

% calculate number of OVERSPILL lakes draining on all days
days_OS_drainages_sum = nansum(days_OS_drainages2,1);
% calculate error in each lake-drainage date (this is the window size)
days_OS_drainages_window = nansum(days_OS_drainages2,2);
average_OS_window = nanmean(days_OS_drainages_window);
% total OVERSPILL drainages in 2022
drainages_OS = length(days_OS_drainages2(:,1))-3 % number of OVERSPILL drainages in 2022 less three STRAGGLERS
drainages_start_OS = 177; 
% drainages_stop_OS = 251;
drainages_stop_OS = 220; 
days_OS = drainages_stop_OS-drainages_start_OS % day window HF drainages occur on in 2022
lambda_OS = drainages_OS./days_OS; 
% lake elevations
elevations_OS = lake_elev(id_OS_drainages);

% log cluster sizes and how many days it takes to make a cluster that large
[~, ~, contributing_rows, cluster_num_lakes,...
    cluster_num_timepoints] = count_simultaneous_drainage(days_OS_drainages2);
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
unique_OS_cluster_S2_2022(:,1) = cluster_time_window_S2_2022(uniqueRowNumbers);
unique_OS_cluster_S2_2022(:,2) = cluster_size_S2_2022(uniqueRowNumbers);

%% %%%%%%%%% v_crit and lake formation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get v_crit and lake-formation time windows
days_OS_form = daily_lake_classifier(id_OS_drainages,:); % HF drainages timeseries
for i=1:length(days_OS_form(:,1))
    % whats the first image we have of this lake?
    day_lake_form_last(i,1) = lake_typeing(id_OS_drainages(i,1),2); % image date of lake first formed
    index1 = find(all_JD_nums(:,1)==day_lake_form_last(i,1)); % Sentinel image index on that day
    if index1>1
    index2 = find(imageID_2022.lake(id_OS_drainages(i),1).plot_IDs>=index1,1,'first'); % Sentinel image on day lake first formed (lake specific)
    % what's the last image available before we see a lake?  (date of last dry basin) 
    index3 = imageID_2022.lake(id_OS_drainages(i),1).plot_IDs(index2-1); % image index before lake first formed (lake specific)
    day_lake_form_first(i,1) = all_JD_nums(index3,1)+1; % date one day after image of dry lake basin (lake specific)
    error_day_lake_form(i,1) = day_lake_form_last(i,1)-day_lake_form_first(i,1); % days of spread (count 24-hr periods)
    else
    day_lake_form_first(i,1) = NaN; % don't plot formation window if it happens before the image timeseries begins
    end
    % v_crit model estimate
    if any(days_OS_form(i,:)==1)
    day_OS_possible(i,1) = find(days_OS_form(i,:)==1,1,'first'); % first day of lake being HFpossible
    else
    day_OS_possible(i,1) = find(days_OS_form(i,:)==2,1,'first'); % first day of lake being HFpossible if classifier goes directly from 0->2
    end
end
% for plotting
formation_average = (day_lake_form_first+day_lake_form_last)/2;
formation_error = day_lake_form_last-formation_average;

%% OS ONLY -- 2022 
figure(2); clf; 
set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',2.0.*[0 0 7.3 13]);
axe1 = axes('Position',[0.10 0.54 0.88 0.44],'Box','on','NextPlot','add');
axe2 = axes('Position',[0.665 0.55 0.30 0.13],'Box','on','NextPlot','add');

fontsize = 10; 
lake_blue = [0    0.4470    0.7410   0.08]; % fourth value is FaceAlpha built in!

axes(axe1); hold on
yyaxis left
i=2; % for legend
errorbar(formation_average, elevations_OS,[],[],formation_error,formation_error,'k.',...
    'LineWidth',0.5,'Marker','none','CapSize',2); % lake formation
plot(day_OS_possible,elevations_OS,'ok',...
    'markerFaceColor','w','MarkerSize',3); % v_crit
plot(1:1:365,days_OS_drainages2(i,:)+elevations_OS(i)-days_OS_drainages2(i,:),...
        'ok','markerFaceColor',lake_blue(1:3),'MarkerSize',3) %[173 221 142]./255);
% RACMO : passing melt thresholds at elevation bands
load('runoff_elev_2022_300m.mat')
melt_threshold = [50, 250, 450, 650, 850, 1050, 1250]; % 50 mm = 5 cm
for i=1:1:length(runoff_elev_2022_300m.elev_bands)-1
    for j=1:length(melt_threshold)
    if runoff_elev_2022_300m.cumulative(i,end) >= melt_threshold(j)
    DOY_threshold(i,j) = find(runoff_elev_2022_300m.cumulative(i,:)>melt_threshold(j),1,'first'); % first day above 0.25 m
    else
    DOY_threshold(i,j) = 1e9;
    end
   end
end
melt_colors = [[248, 250, 165],;
       [248, 239, 144],;
       [249, 228, 124],;
       [249, 217, 104],;
       [250, 206, 83],;
       [250, 195, 63],;
       [251, 185, 43]]./255;
area(DOY_threshold(:,1),runoff_elev_2022_300m.elev_bands(1:end-1)+25,'EdgeColor','none','FaceColor',melt_colors(1,:)); 
area(DOY_threshold(:,2),runoff_elev_2022_300m.elev_bands(1:end-1)+25,'EdgeColor','none','FaceColor',melt_colors(2,:)); 
area(DOY_threshold(:,3),runoff_elev_2022_300m.elev_bands(1:end-1)+25,'EdgeColor','none','FaceColor',melt_colors(3,:)); 
area(DOY_threshold(:,4),runoff_elev_2022_300m.elev_bands(1:end-1)+25,'EdgeColor','none','FaceColor',melt_colors(4,:)); 
area(DOY_threshold(:,5),runoff_elev_2022_300m.elev_bands(1:end-1)+25,'EdgeColor','none','FaceColor',melt_colors(5,:)); 
area(DOY_threshold(:,6),runoff_elev_2022_300m.elev_bands(1:end-1)+25,'EdgeColor','none','FaceColor',melt_colors(6,:));
area(DOY_threshold(:,7),runoff_elev_2022_300m.elev_bands(1:end-1)+25,'EdgeColor','none','FaceColor',melt_colors(7,:));
% labels
text(269,1650,{'Runoff:'},'FontSize',fontsize+1,'FontName','Avenir','Rotation',0,'HorizontalAlignment','right');
text(269,1560,sprintf('%d cm',melt_threshold(1)./10),'FontName','Avenir','Rotation',0,'HorizontalAlignment','right');
text(269,1470,sprintf('%d cm',melt_threshold(2)./10),'FontName','Avenir','Rotation',0,'HorizontalAlignment','right');
text(269,1400,sprintf('%d cm',melt_threshold(3)./10),'FontName','Avenir','Rotation',0,'HorizontalAlignment','right');
text(269,1350,sprintf('%d cm',melt_threshold(4)./10),'FontName','Avenir','Rotation',0,'HorizontalAlignment','right');
text(269,1290,sprintf('%d cm',melt_threshold(5)./10),'FontName','Avenir','Rotation',0,'HorizontalAlignment','right');
text(269,1210,sprintf('%d cm',melt_threshold(6)./10),'FontName','Avenir','Rotation',0,'HorizontalAlignment','right');
text(269,1140,sprintf('%d cm',melt_threshold(7)./10),'FontName','Avenir','Rotation',0,'HorizontalAlignment','right');

errorbar(formation_average, elevations_OS,[],[],formation_error,formation_error,'k.',...
    'LineWidth',0.5,'Marker','none','CapSize',2); % lake formation
plot(day_OS_possible,elevations_OS,'ok',...
    'markerFaceColor','w','MarkerSize',3); % v_crit

% OS drainages
for i=1:1:length(days_OS_drainages2(:,1))
    plot(1:1:365,days_OS_drainages2(i,:)+elevations_OS(i)-days_OS_drainages2(i,:),'ok',...
        'markerFaceColor',lake_blue(1:3),'MarkerSize',3); 
end
xlim([154 270]); ylim([600 1800]);
set(gca,'ytick',500:100:1800);
set(gca,'LineWidth',1.1,'tickdir','in','FontSize',fontsize,'FontName','Avenir');
set(gca,'ycolor','k');
%ylabel('Ice-sheet Surface Elevation [ m a.s.l. ]','FontSize',fontsize);
xlabel('Day of Year, 2022','FontSize',fontsize)
%title('Overspill events constrained by Sentinel-2','FontSize',fontsize)
set(gca, 'SortMethod', 'depth')
text(141.4, 1800, 'e','FontName','Arial','FontWeight','bold','FontSize',fontsize+4);
lgd=legend('Formation','{\it v_{crit}}','Overspill event','Location','NorthEast','Orientation','vertical');
set(lgd,'Box','on','FontName','Avenir','FontSize',fontsize-1,'LineWidth',0.5,...
    'Position',[0.224 0.941 0.01 0.01])

yyaxis right
hold on;
xlim([154 270]); ylim([0 max(days_OS_drainages_sum)+1]); set(gca,'ycolor','k');
set(gca, 'SortMethod', 'depth','YTickLabel',[],'YTick',[])
grid on; box on
ax=gca; ax.Layer = 'top';
for i=1:1:length(m)
    for j=1:1:length(t)
    PDF(i,j) = (((lambda_OS.*t(j)).^(m(i))).*(exp(-1.*lambda_OS.*t(j))))./...
        (factorial(m(i)));       
    end
end
% probability of observing one or more events less than or equal to size m = 
% 1 minus the probability of observing zero events + an event of size 1 +
% + an event of size 2 + ... + an event of size (m-1) + an event of size m:
CDF_2 = 1-cumsum(PDF,1);

% find values at a given \alpha threshold: (~1.23 * t * \lambda)
for i=1:1:length(t)
    [sig_OS_cluster(i,1)] = find(CDF_2(:,i)<=alpha,1)-1;
end

% clusters observed in 2022, as a f(cluster size, cluster probability),
% given the time window observed
for i=1:1:length(unique_OS_cluster_S2_2022(:,1))
    cluster_probability_OS_S2_2022(i,1) = CDF_2(unique_OS_cluster_S2_2022(i,2)+1,... % num of events
                                unique_OS_cluster_S2_2022(i,1)); % time window
end
unique_OS_cluster_S2_2022(unique_OS_cluster_S2_2022(:,1)>=40,1:2) = NaN;

axes(axe2); 
hold on;
plot(0:1:100,(0:1:100).*lambda_OS,'-','Color',[0.8500    0.3250    0.0980],'LineWidth',1.2);
a1 = area(0:1:length(sig_OS_cluster),vertcat(0,sig_OS_cluster));
a1.FaceColor = [0.8500    0.3250    0.0980]; a1.FaceAlpha = 0.25;
a1.FaceColor = lake_blue(1:3); a1.FaceAlpha = 0.15;
a1.EdgeColor = [0.8500    0.3250    0.0980]; a1.BaseLine.LineStyle = 'none';
a1.BaseValue = 100; a1.LineWidth = 1.9; a1.LineStyle = ':';
plot(unique_OS_cluster_S2_2022(:,1),unique_OS_cluster_S2_2022(:,2),'ok',... 
    'MarkerFaceColor','none','MarkerSize',4); 
% significant clusters
for i=1:length(unique_OS_cluster_S2_2022(:,1))
    if cluster_probability_OS_S2_2022(i) <= alpha
        plot(unique_OS_cluster_S2_2022(i,1),unique_OS_cluster_S2_2022(i,2),'kp',...
        'MarkerFaceColor',[247 252 185]./255,'MarkerSize',9);
    else
    end
end
text(1,52,'f','FontName','Arial','FontWeight','bold','FontSize',14)
n_unique_clusters_2022 = length(unique_OS_cluster_S2_2022(:,1))
text(14.5,40,sprintf('\\lambda = %3.2f',lambda_OS),'FontName','Arial','FontSize',10,'rotation',38,...
    'Color',[0.8500    0.3250    0.0980],'FontWeight','bold');
ylabel('No. of OS events, {\itm}'); xlabel('Time window, {\itw} [ days ]'); 
ylim([0 60]); xlim([0 20]); grid on; box on;
set(gca,'LineWidth',1.1,'tickdir','in','FontSize',9,'XAxisLocation','top');

%% 2023 
clear all;

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

% ROI that's relevant to the GNSS array, after FASTER locations update (2025-April):
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

% polygon ROI check
%figure; plot(ROI_relevant_xv,ROI_relevant_yv)
%hold on; plot(environs_lakes.X_km_origin,environs_lakes.Y_km_origin,'o')
%title('2023'); %pause; close all

% racmo at time of HF drainage
load cumulative_runoff_last_day_to_HF_2023.mat

% preliminaries for homo poisson
alpha = 0.05; % significance threshold
t = 1:1:60; % time window [ days ]
m = 0:1:250; % number of events [ cluster size ]

% drainage dates for each drainage type OVERSPILL
id_OS_drainages = find(lake_typeing(:,7)==3); % OVERSPILL drainages row IDs [ OS == 3 ]
id_OS_drainages = id_OS_drainages(1:end-1); % 
days_OS_drainages = daily_lake_classifier(id_OS_drainages,:); % OS drainages timeseries
days_OS_drainages(days_OS_drainages~=2) = NaN; % NaN if not in filling or drainage window
% remove rows that are all NaNs
days_OS_drainages = days_OS_drainages(~all(isnan(days_OS_drainages), 2), :);
%figure(3); imagesc(days_OS_drainages); colorbar; xlim([160 260]); shg
days_OS_drainages = days_OS_drainages-1;

% % get the in-between days, and not the initial day, --> 1
days_OS_drainages2 = NaN(length(id_OS_drainages),365); % prefill with NaNs
for i=1:1:length(id_OS_drainages)
    % if just one day draining window: it's that day
    if environs_lakes.laketypeing_dates(id_OS_drainages(i),3) == environs_lakes.laketypeing_dates(id_OS_drainages(i),4)
        days_OS_drainages2(i,environs_lakes.laketypeing_dates(id_OS_drainages(i),4)) = 1;
    elseif environs_lakes.laketypeing_dates(id_OS_drainages(i),3)+1 == environs_lakes.laketypeing_dates(id_OS_drainages(i),4)
        days_OS_drainages2(i,environs_lakes.laketypeing_dates(id_OS_drainages(i),4)) = 1;
    else 
        days_OS_drainages2(i,(environs_lakes.laketypeing_dates(id_OS_drainages(i),3)+1):environs_lakes.laketypeing_dates(id_OS_drainages(i),4)) = 1;
   end
end

%figure(4); imagesc(days_OS_drainages2); colorbar; xlim([160 260]); shg
% calculate number of OVERSPILL lakes draining on all days
days_OS_drainages_sum = nansum(days_OS_drainages2,1);
% calculate error in each lake-drainage date (this is the window size)
days_OS_drainages_window = nansum(days_OS_drainages2,2);
average_OS_window = nanmean(days_OS_drainages_window);
% total OVERSPILL drainages in 2022
drainages_OS = length(days_OS_drainages2(:,1)) % number of OVERSPILL drainages in 2023
drainages_start_OS = 177; 
% drainages_stop_OS = 240;
drainages_stop_OS = 215;
days_OS = drainages_stop_OS-drainages_start_OS % day window HF drainages occur on in 2023
lambda_OS = drainages_OS./days_OS; 
% lake elevations
elevations_OS = lake_elev(id_OS_drainages);

% log cluster sizes and how many days it takes to make a cluster that large
[num_drainage_timepoints, drainage_indices, contributing_rows, cluster_num_lakes,...
    cluster_num_timepoints] = count_simultaneous_drainage(days_OS_drainages2);
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
unique_OS_cluster_S2_2022(:,1) = cluster_time_window_S2_2022(uniqueRowNumbers);
unique_OS_cluster_S2_2022(:,2) = cluster_size_S2_2022(uniqueRowNumbers);

%% %%%%%%%%% v_crit and lake formation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get v_crit and lake-formation time windows
days_OS_form = daily_lake_classifier(id_OS_drainages,:); % HF drainages timeseries
for i=1:length(days_OS_form(:,1))
    % whats the first image we have of this lake?
    day_lake_form_last(i,1) = lake_typeing(id_OS_drainages(i,1),2); % image date of lake first formed
    index1 = find(all_JD_nums(:,1)==day_lake_form_last(i,1)); % Sentinel image index on that day
    if index1>1
    index2 = find(imageID_2023.lake(id_OS_drainages(i),1).plot_IDs>=index1,1,'first'); % Sentinel image on day lake first formed (lake specific)
    % what's the last image available before we see a lake?  (date of last dry basin) 
        if index2>1
        index3 = imageID_2023.lake(id_OS_drainages(i),1).plot_IDs(index2-1); % image index before lake first formed (lake specific)
        day_lake_form_first(i,1) = all_JD_nums(index3,1)+1; % date one day after image of dry lake basin (lake specific)
        else
        day_lake_form_first(i,1) = all_JD_nums(index2,1)-1; % plot 1 days back if no prior image
        end
    error_day_lake_form(i,1) = day_lake_form_last(i,1)-day_lake_form_first(i,1); % days of spread (count 24-hr periods)
    else
    day_lake_form_first(i,1) = NaN; % don't plot formation window if it happens before the image timeseries begins
    end
    % v_crit model estimate
    if any(days_OS_form(i,:)==1)
    day_OS_possible(i,1) = find(days_OS_form(i,:)==1,1,'first'); % first day of lake being HFpossible
    else
    day_OS_possible(i,1) = find(days_OS_form(i,:)==2,1,'first'); % first day of lake being HFpossible if classifier goes directly from 0->2
    end
end
% for plotting
formation_average = (day_lake_form_first+day_lake_form_last)/2;
formation_error = day_lake_form_last-formation_average;

%% plotting 
axe3 = axes('Position',[0.10 0.045 0.88 0.44],'Box','on','NextPlot','add');
axe4 = axes('Position',[0.665 0.29 0.30 0.183],'Box','on','NextPlot','add');
fontsize = 10; 
lake_blue = [0    0.4470    0.7410   0.08]; % fourth value is FaceAlpha built in!

axes(axe3); hold on
yyaxis left
% RACMO : passing melt thresholds at elevation bands
load('runoff_elev_2023_300m.mat')
melt_threshold = [50, 250, 450, 650, 850, 1050, 1250, 1450, 1650, 1850]; 
for i=1:1:length(runoff_elev_2023_300m.elev_bands)-1
    for j=1:length(melt_threshold)
    if runoff_elev_2023_300m.cumulative(i,end) >= melt_threshold(j)
    DOY_threshold(i,j) = find(runoff_elev_2023_300m.cumulative(i,:)>melt_threshold(j),1,'first'); % first day above 0.25 m
    else
    DOY_threshold(i,j) = 1e9;
    end
   end
end
melt_colors = [[248, 250, 165],;
       [248, 239, 144],;
       [249, 228, 124],;
       [249, 217, 104],;
       [250, 206, 83],;
       [250, 195, 63],;
       [251, 185, 43]]./255;
area(DOY_threshold(:,1),runoff_elev_2023_300m.elev_bands(1:end-1),'EdgeColor','none','FaceColor',melt_colors(1,:)); 
area(DOY_threshold(:,2),runoff_elev_2023_300m.elev_bands(1:end-1),'EdgeColor','none','FaceColor',melt_colors(2,:)); 
area(DOY_threshold(:,3),runoff_elev_2023_300m.elev_bands(1:end-1),'EdgeColor','none','FaceColor',melt_colors(3,:)); 
area(DOY_threshold(:,4),runoff_elev_2023_300m.elev_bands(1:end-1),'EdgeColor','none','FaceColor',melt_colors(4,:)); 
area(DOY_threshold(:,5),runoff_elev_2023_300m.elev_bands(1:end-1),'EdgeColor','none','FaceColor',melt_colors(5,:)); 
area(DOY_threshold(:,6),runoff_elev_2023_300m.elev_bands(1:end-1),'EdgeColor','none','FaceColor',melt_colors(6,:));
area(DOY_threshold(:,7),runoff_elev_2023_300m.elev_bands(1:end-1),'EdgeColor','none','FaceColor',melt_colors(7,:));
% vcrit, formation dates
errorbar(formation_average, elevations_OS,[],[],formation_error,formation_error,'k.',...
    'LineWidth',0.5,'Marker','none','CapSize',2); % lake formation
plot(day_OS_possible,elevations_OS,'ok',...
    'markerFaceColor','w','MarkerSize',3); % v_crit

% OS drainages
for i=1:1:length(days_OS_drainages2(:,1))
    plot(1:1:365,days_OS_drainages2(i,:)+elevations_OS(i)-days_OS_drainages2(i,:),'ok',...
        'markerFaceColor',lake_blue(1:3),'MarkerSize',3);
end
xlim([154 270]); ylim([600 1800]);
set(gca,'ytick',600:100:1800);
set(gca,'LineWidth',1.1,'tickdir','in','FontSize',fontsize,'FontName','Avenir');
set(gca,'ycolor','k');
%ylabel('Ice-sheet Surface Elevation [ m a.s.l. ]','FontSize',fontsize);
xlabel('Day of Year, 2023','FontSize',fontsize)
set(gca, 'SortMethod', 'depth')
text(141.4, 1800, 'g','FontName','Arial','FontWeight','bold','FontSize',fontsize+4);

yyaxis right
hold on;
xlim([154 270]); ylim([0 max(days_OS_drainages_sum)]);
set(gca,'ycolor','k');
set(gca, 'SortMethod', 'depth','YTickLabel',[],'YTick',[])
grid on; box on
ax=gca; ax.Layer = 'top';
for i=1:1:length(m)
    for j=1:1:length(t)
    PDF(i,j) = (((lambda_OS.*t(j)).^(m(i))).*(exp(-1.*lambda_OS.*t(j))))./...
        (factorial(m(i)));

    end
end
% probability of observing one or more events less than or equal to size m = 
% 1 minus the probability of observing zero events + an event of size 1 +
% + an event of size 2 + ... + an event of size (m-1) + an event of size m:
CDF_2 = 1-cumsum(PDF,1);
% find values at a given \alpha threshold: (~1.23 * t * \lambda)
for i=1:1:length(t)
    [sig_OS_cluster(i,1)] = find(CDF_2(:,i)<=alpha,1)-1;
end
% clusters observed in 2022, as a f(cluster size, cluster probability),
% given the time window observed
for i=1:1:length(unique_OS_cluster_S2_2022(:,1))
    cluster_probability_OS_S2_2022(i,1) = CDF_2(unique_OS_cluster_S2_2022(i,2)+1,... % num of events
                                unique_OS_cluster_S2_2022(i,1)); % time window
end

axes(axe4); 
hold on;
plot(0:1:100,(0:1:100).*lambda_OS,'-','Color',[0.8500    0.3250    0.0980],'LineWidth',1.2);
a1 = area(0:1:length(sig_OS_cluster),vertcat(0,sig_OS_cluster));
a1.FaceColor = [0.8500    0.3250    0.0980]; a1.FaceAlpha = 0.25;
a1.FaceColor = lake_blue(1:3); a1.FaceAlpha = 0.15;
a1.EdgeColor = [0.8500    0.3250    0.0980]; a1.BaseLine.LineStyle = 'none';
a1.BaseValue = 150; a1.LineWidth = 1.9; a1.LineStyle = ':';
plot(unique_OS_cluster_S2_2022(:,1),unique_OS_cluster_S2_2022(:,2),'ok',... 
    'MarkerFaceColor','none','MarkerSize',4); 
% significant clusters
for i=1:length(unique_OS_cluster_S2_2022(:,1))
    if cluster_probability_OS_S2_2022(i) <= alpha
        plot(unique_OS_cluster_S2_2022(i,1),unique_OS_cluster_S2_2022(i,2),'kp',...
        'MarkerFaceColor',[247 252 185]./255,'MarkerSize',9);
    else
    end
end
text(1,72,'h','FontName','Arial','FontWeight','bold','FontSize',14)
n_unique_clusters_2023 = length(unique_OS_cluster_S2_2022(:,1))
text(14,43,sprintf('\\lambda = %3.2f',lambda_OS),'FontName','Arial','FontSize',10,'rotation',40,...
    'Color',[0.8500    0.3250    0.0980],'FontWeight','bold');
ylabel('Number of OS events, {\itm}'); xlabel('Time window, {\itw} [ days ]'); 
ylim([0 80]); xlim([0 20]); grid on; box on;
set(gca,'LineWidth',1.1,'tickdir','in','FontSize',9);


%% print figure
print(gcf,'-dpng','-r600',sprintf('../paperfigs/paperfig4_homopo_vcrit_OS_260128_repository.png')); 
%exportgraphics(gcf,'../paperfigs/paperfigED4_homopo_vcrit_OS_250528_repository.pdf',"ContentType","vector");