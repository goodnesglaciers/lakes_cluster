function [band_var] = nevis_elevband(fn,tis) % [uncomment to make a function ]
% nevis_elevband(fn,tis)
% 
% plots discharge, height of cavity layer, effective pressure, and area of
% channels through time as a function of 100-m ice surface elevation bins. 
% 
% 24 January 2017: created for a way to compare nevis output to SAR 
% velocity maps (LAS) 

%% load initial timestep
if exist([fn,'/0000']), load([fn,'/0000']); else load(fn); end
load([fn,'/0001']);
if isfield(oo,'reversey'), reversey = oo.reversey; else reversey = 0; end

%% contour levels
% z_conts = -1000:50:5000;
z_conts = -1000:10:5000;
b_conts = z_conts;
s_conts = z_conts;
p_conts = z_conts/10;

%% axes limits
axx = (ps.x/10^3)*[gg.xl-0.01 gg.xr gg.yb-0.01 gg.yt];

%% extract initial variables
if isfield(vv,'nbdy'), gg = nevis_label(gg,vv.nbdy); aa.phi = aa.phi_a(gg.nbdy); end % re-mask
oo.evaluate_variables = 1; oo.evaluate_residual = 0; oo.evaluate_jacobian = 0; [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
nevis_unpack(aa,gg);

%get rid of points outside domain
nx(gg.nout) = NaN;
ex(gg.eout) = NaN;
fx(gg.fout) = NaN;
cx(gg.cout) = NaN;
    
%boundary curve
x_out = gg.nx(gg.n1); y_out = gg.ny(gg.n1);
if ~isempty(x_out), tmp = nevis_orderboundary(x_out,y_out); x_out = x_out(tmp); y_out = y_out(tmp); end % reorder to follow boundary

%assign surface elevation into bins
s_bands = (0:100:1600)./ps.z;           % elevation bins
s_disc = discretize2(aa.s,s_bands);     % discretize surface elvation into bins
s_disc = reshape(s_disc,gg.nI,gg.nJ);

%% load timesteps
for i_t = 1:length(tis)
disp(['Time step ',num2str(i_t),' / ',num2str(length(tis)),' ...'])

%% load timestep
load([fn,'/',int2four(tis(i_t))]);

%% extract new variables
if isfield(vv,'nbdy'), gg = nevis_label(gg,vv.nbdy); aa.phi = aa.phi_a(gg.nbdy); end % remask
aa = nevis_inputs(vv.t,aa,pp,gg,oo);
oo.evaluate_variables = 1; oo.evaluate_residual = 0; oo.evaluate_jacobian = 0; [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
nevis_unpack(aa,vv2);
clear tt vv vv2

% extract variables of interest
    % discharge
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny; 
    zz_q = ps.qs*reshape(qs+qe+qQ,gg.nI,gg.nJ);
    cax = [10^(-5) 10^(0)]; 
    
    % effective pressure
    zz_N = (ps.phi/10^6)*reshape(phi_0-phi,gg.nI,gg.nJ); 
    cax = [-1 1]; % [-3 3]; 
    
    % elastic sheet 
    zz_he = (ps.h*100)*reshape(he,gg.nI,gg.nJ); 
    cax = [0 100];
    
    % sheet thickness
    zz_hs = (ps.h*100)*reshape(hs,gg.nI,gg.nJ); 
    cax = [0 50];
    
    % surface input
    zz_E = (ps.m*10^3*pd.td)*reshape(E,gg.nI,gg.nJ); 
    cax = [0 80]; 
    
% discretize variables    
for j=1:(length(s_bands)-1)
    zzq_disc2 = s_disc;  zzq_disc2(zzq_disc2 ~= j)= NaN;    % NaN except for the band of interest
    zzq_disc3 = zzq_disc2./zzq_disc2;                       % 1 (of interest); NaN (not of interest)
    
    interest_zzq = (zzq_disc3).*zz_q;
    band_q(i_t,j) = nanmean(nanmean(interest_zzq));
    
    interest_zzN = (zzq_disc3).*zz_N;
    band_N(i_t,j) = nanmean(nanmean(interest_zzN));
    
    interest_zzhe = (zzq_disc3).*zz_he;
    band_he(i_t,j) = nanmean(nanmean(interest_zzhe));
    
    interest_zzhs = (zzq_disc3).*zz_hs;
    band_hs(i_t,j) = nanmean(nanmean(interest_zzhs));   
    
    interest_zzE = (zzq_disc3).*zz_E;
    band_E(i_t,j) = nanmean(nanmean(interest_zzE));
end

% Variables into structures
band_var.band_q = band_q;   band_var.band_N = band_N;
band_var.band_he = band_he; band_var.band_hs = band_hs;
band_var.band_E = band_E;

end

%% FIGURE OVER TIME

figure(10)
clf

t = tis;
elev = (s_bands(1:16).*ps.z) + 50;
c_time = parula(365);

for i=1:length(t)
    ax(3) = subplot(3,1,3);        
        plot(elev,-1.*band_N(i,:),'Color',c_time(t(i),:));
        hold on;
        xlabel('Elevation [ m ]'); xlim([0 1600]);
        ylabel('-N [ MPa ]');
        ylim([min(min(-1.*band_N)) max(max(-1.*band_N))])
        
    ax(2) = subplot(3,1,2); 
        %for i=1:length(t)
        plot(elev,band_q(i,:),'Color',c_time(t(i),:));
        hold on;              
        xlabel('Elevation [ m ]'); xlim([0 1600]);
        ylabel('Q [ m^3/s ]');
        ylim([ min(min(band_q)) max(max(band_q)) ])
        
    ax(1) = subplot(3,1,1); 
        %for i=1:length(t)
        plot(elev,band_E(i,:),'Color',c_time(t(i),:));
        hold on;
        xlabel('Elevation [ m ]'); xlim([0 1600]);
        ylabel('E [ m^3/s ]');        
end

% color bar
   %cbar_size = [0.4 0.9 0.2 0.02];
   cx = colorbar('peer',gca,'horizontal','position',[0.4 0.95 0.2 0.01],'XAxislocation','top');
   %cx.Label.String = 't [ d ]'; cx.Label.Units = 'normalized'; cx.Label.Position = [-0.05 0]; cx.Label.VerticalAlignment = 'bottom'; cx.Label.HorizontalAlignment = 'right';
   colormap(parula(365));
   
   %text(10,10,'t [ d ]','FontSize',12,'Interpreter','tex');





