function [run] = runoff_function(t, gg)
% 10 May 2016 
% Laura Stevens
% load RACMO daily runoff values over 160 km region

load('R1km_for_nevis.mat')         % load data (previously collated)
day = round(t);                    % round time to nearest interger day
year = sprintf('y%d',gg.year);     % year of interest in string
run = R1km.(year).run(:,:,day);    % runoff on year and day of interest
end