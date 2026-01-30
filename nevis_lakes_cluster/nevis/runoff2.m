function [r] = runoff(day,runoff_2011_nevis)
% 10 May 2016 LAS
% load RACMO daily runoff values over 160 km region
day = round(day);                    % round time to nearest interger day
% reshape to gg.nIJ (distributed input at nodes)
r = ((runoff_2011_nevis(day,:))./(1000*86400))'; % runoff on year and day of interest (m/sec)
end