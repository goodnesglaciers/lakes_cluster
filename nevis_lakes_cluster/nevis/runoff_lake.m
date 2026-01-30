function [r] = runoff_lake(day, runoff_2011_nevis, runoff_time)
% % 10 May 2016 LAS
% day = round(day);       % round time to nearest interger day

% % 19 Dec 2019 LAS
% find index of closest time in runoff_time to model step day
[~,idx] = min(abs(day-runoff_time));

% mm w.e./day  to m/s
r = ((runoff_2011_nevis(idx,:))/(1000*86400))'; % runoff on day of interest (m/sec)
end