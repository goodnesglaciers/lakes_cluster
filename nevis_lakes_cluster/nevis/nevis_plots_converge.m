%% nevis plots for converge shape
%% 11 November 2015 LAS
time = [0:.1:40];

%% set up figure
fig1=figure(1); clf; set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',[5 2 15 12]);
set(0,'DefaultAxesFontSize',10);
set(0,'DefaultTextFontSize',10);
axes_size = [0.1 0.1 0.8 0.8];
axe1 = axes('Position',[ 0.12 0.1 0.8 0.8 ],'Box','on','NextPlot','add');

axes(axe1)
% 51 nodes
load nevis_151111a.mat
for i=1:length(time)
[row(i),col(i)] = find(tt(1,i).pts_N < 0, 1, 'last');
end
plot(time,gg.nx(row,1)*ps.x/10e3,'k') 
hold all

load nevis_151111b.mat
for i=1:length(time)
[row(i),col(i)] = find(tt(1,i).pts_N < 0, 1, 'last');
end
plot(time,gg.nx(row,1)*ps.x/10e3,'b') 
hold all

load nevis_151111c.mat
for i=1:length(time)
[row(i),col(i)] = find(tt(1,i).pts_N < 0, 1, 'last');
end
plot(time,gg.nx(row,1)*ps.x/10e3,'g') 
hold all

load nevis_151111d.mat
for i=1:length(time)
[row(i),col(i)] = find(tt(1,i).pts_N < 0, 1, 'last');
end
plot(time,gg.nx(row,1)*ps.x/10e3,'r') 
hold all

xlabel(' Time [ days ] '); ylabel(' X-location [ km ] ');
title(' Position of last value where N < 0');
legend('51 nodes','101 nodes','151 nodes','201 nodes');
