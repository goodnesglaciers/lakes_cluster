%function nevis_plot9(vv,aa,pp,ps,gg,oo) % [uncomment to make a function ]
% nevis_plot5(vv,aa,pp,ps,gg,oo)
% comment first line to use as script taking inputs from current workspace
% plots discharge, height of cavity layer, effective pressure, and area of
% channels, topography, ice thickness, surface velocity, and englacial 
% storage, and supraglacial inputs in a 9-panel plot. 
% 
% 28 July 2016: taken from nevis_channelx for regional channel cross sections (LAS) 
% 28 March 2019: taken from nevis_plot4 for Helheim model (LAS)

%% options
fn = [oo.root,oo.fn];
if isfield(oo,'save_plot'), save_plot = oo.save_plot; else save_plot = 0; end
if isfield(oo,'reversey'), reversey = oo.reversey; else reversey = 0; end

%% figures to plot
if isfield(oo,'discharge'), discharge = oo.discharge; else discharge = 1; end
if isfield(oo,'topography'), topography = oo.topography; else topography = 0; end
if isfield(oo,'discharge_lines'), discharge_lines = oo.discharge_lines; else discharge_lines = 0; end
if isfield(oo,'velocity'), velocity = oo.velocity; else velocity = 0; end
if isfield(oo,'input'), input = oo.input; else input = 0; end
if isfield(oo,'dissipation'), dissipation = oo.dissipation; else dissipation = 0; end
if isfield(oo,'thickness'), thickness = oo.thickness; else thickness = 0; end
if isfield(oo,'sheet'), sheet = oo.sheet; else sheet = 0; end
if isfield(oo,'elastic'), elastic = oo.elastic; else elastic = 0; end
if isfield(oo,'area'), area = oo.area; else area = 0; end
if isfield(oo,'pressure'), pressure = oo.pressure; else pressure = 0; end
if isfield(oo,'channelx'), channelx = oo.channelx; else channelx = 0; end

%% elevation range
z_range = [-999 1000];
% z_range = [2000 3000];

%% contour levels
z_conts = -1000:100:5000;
b_conts = z_conts;
s_conts = z_conts;
p_conts = z_conts/100;

%% extract variables
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
if ~isempty(gg.n1);
x_out = gg.nx(gg.n1); y_out = gg.ny(gg.n1);
tmp = nevis_orderboundary(x_out,y_out); x_out = x_out(tmp); y_out = y_out(tmp); % reorder to follow boundary
else x_out = []; y_out = [];
end

%% axes limits
axx = (ps.x/10^3)*[gg.xl gg.xr gg.yb gg.yt]; axx(3) = -5.1;

%% set up figure
figure(1); set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',[1 1 170 130]./5);
set(0,'DefaultAxesFontSize',10);
set(0,'DefaultTextFontSize',10);
axes_size = [0.1 0.1 0.4 0.4];

axe1 = axes('Position',[0.02 0.67 0.3 0.3],'Box','on','NextPlot','add','XTickLabel',[]); 
axe2 = axes('Position',[0.02 0.36 0.3 0.3],'Box','on','NextPlot','add','XTickLabel',[]); 
axe3 = axes('Position',[0.02 0.05 0.3 0.3],'Box','on','NextPlot','add'); 
axe4 = axes('Position',[0.35 0.67 0.3 0.3],'Box','on','NextPlot','add','XTickLabel',[],'YTickLabel',[]); 
axe5 = axes('Position',[0.35 0.36 0.3 0.3],'Box','on','NextPlot','add','XTickLabel',[],'YTickLabel',[]); 
axe6 = axes('Position',[0.35 0.05 0.3 0.3],'Box','on','NextPlot','add','YTickLabel',[]); 
axe7 = axes('Position',[0.67 0.67 0.3 0.3],'Box','on','NextPlot','add','XTickLabel',[],'YTickLabel',[]); 
axe8 = axes('Position',[0.67 0.36 0.3 0.3],'Box','on','NextPlot','add','XTickLabel',[],'YTickLabel',[]); 
axe9 = axes('Position',[0.67 0.05 0.3 0.3],'Box','on','NextPlot','add','YTickLabel',[]); 

%%
%adjust axes positions to account for aspect ratio of figure
aspect = (axx(4)-axx(3))/(axx(2)-axx(1)); % aspect of axes limits
pos = axes_size; tmp = get(gcf,'position'); sca = tmp(4)/tmp(3);
if aspect>(pos(4)/pos(3))*sca, 
    pos = pos + [pos(4)*sca*(1-1/aspect)/2 0 -pos(3)+pos(4)*sca*(1/aspect) 0]; % constrained by y
else
    pos = pos + [0 pos(3)/sca*(1-aspect)/2 0 -pos(4)+pos(3)/sca*aspect]; % constrained by x
end 
axes_size = pos;

% %% load Joughin 2013 The Cryosphere TerraSAR-X footprint
% radius=6378137.0; eccen=0.08181919; lat_true=70; lon_posy=-45; % projection parameters
% foot = [68.884795, -50.182106; 68.952545, -49.409098;
%         %68.393966, -48.867272; 68.282561, -49.831050];
%         68.469744, -48.970285; 68.399076, -49.709055];
% moulin = [68.723585, -49.536195]; % moulin is (0,0)    
% [moulin_x,moulin_y] = polarstereo_fwd(moulin(1),moulin(2),radius,eccen,lat_true,lon_posy);
% [foot_x,foot_y] = polarstereo_fwd(foot(:,1),foot(:,2),radius,eccen,lat_true,lon_posy);
% footprint(:,1) = foot_x-moulin_x; footprint(:,2) = foot_y - moulin_y;
% 
% FLXX = [68.7735, 310.0752-360;
%         68.7559, 310.2587-360;
%         68.7373, 310.4469-360;
%         68.7253, 310.5542-360;
%         68.7158, 310.8561-360;
%         68.7192, 311.1320-360];
% [flxx_x,flxx_y] = polarstereo_fwd(FLXX(:,1),FLXX(:,2),radius,eccen,lat_true,lon_posy);
% flxx_plot(:,1) = flxx_x-moulin_x; flxx_plot(:,2) = flxx_y - moulin_y; 

%% plot

axes(axe1);
    % bed
    xx = (ps.x/10^3)*gg.nx; yy = (ps.x/10^3)*gg.ny;  % keep outside points here to prevent edge effects disappearing
    zz = ps.z*reshape(b,gg.nI,gg.nJ); zz(zz<=-1500) = NaN;
    cax = [-1500 1500]; 
    
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis     
%     hand = pcolor(xx,yy,zz); set(hand,'linestyle','none'); shading interp;
    hand = imagesc(xx(:,1),yy(1,:),zz'); set(gca,'YDir','normal');
    hold on;
    caxis(cax); 
    axis image;
    axis(axx);
    
    % color bar
    tmp = get(gca,'position'); 
    cbar_size = [tmp(1)+tmp(3)-0.02 tmp(2)+tmp(4)*1/8 0.01 tmp(4)*6/8];
    cx1 = colorbar('peer',gca,'vertical','position',cbar_size,'YAxislocation','right');
    text(5,22,'z [ m ]','VerticalAlignment','bottom','HorizontalAlignment','right');

%     % add bed contours
%     zz = ps.z*reshape(b,gg.nI,gg.nJ); zz(zz<=-400) = NaN;
%     fil = fspecial('average'); zzz = conv2(zz,fil,'valid'); %zzz = zz;
%     xrg = (size(zz,1)-size(zzz,1))/2+(1:size(zzz,1)); yrg = (size(zz,2)-size(zzz,2))/2+(1:size(zzz,2));
%    hold on; contour(xx(xrg,yrg),yy(xrg,yrg),zzz,b_conts,'color',0.5*[1 1 1]);
 
%     % add moulins
%     if ~isfield(pp,'ni_m'), pp.ni_m = []; end
%     mscale = 100;  %amount to scale moulin input by to make size of moulin dot
%     x = (ps.x/10^3)*nx(pp.ni_m);
%     y = (ps.x/10^3)*ny(pp.ni_m);
%     hold on;
%     for i_m = 1:length(pp.ni_m),
%         if E(pp.ni_m(i_m))>0,
%             plot(x(i_m),y(i_m),'ko','Markersize',4+E(pp.ni_m(i_m))/mscale,'MarkerFaceColor',1*[1 1 1]); % mark moulins   
%         else
%             plot(x(i_m),y(i_m),'ko','Markersize',4,'MarkerFaceColor',0.8*[1 1 1]); % mark moulins  
%         end
%     end
    
    % add boundary
    plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5); 
    
    cmap1=load('cmapland2.mat'); colormap(gca,cmap1.cmap); 
    text(-0.11,0.5,'y [ km ]','HorizontalAlignment','center','units','normalized','rotation',90);
    set(gca,'TickDir','out','LineWidth',1.1); box on;
    text(-34,26.5,([oo.root,oo.fn]),'FontSize',12,'Interpreter','latex');
    
axes(axe2);
    % ice thickness
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny;  
    zz = ps.z*reshape(H,gg.nI,gg.nJ); 
    cax = [0 2000]; 
    
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis     
    hand = pcolor(xx,yy,zz); set(hand,'linestyle','none'); % shading interp;
    caxis(cax); 
    axis image;
    axis(axx);
    
    % add surface contours
    zz = ps.z*reshape(s,gg.nI,gg.nJ);
    fil = fspecial('average'); zzz = conv2(zz,fil,'valid'); %zzz = zz;
    xrg = (size(zz,1)-size(zzz,1))/2+(1:size(zzz,1)); yrg = (size(zz,2)-size(zzz,2))/2+(1:size(zzz,2));
    hold on; [CC,hh]=contour(xx(xrg,yrg),yy(xrg,yrg),zzz,s_conts,'color',0.5*[1 1 1],'linewidth',1);
    vv=[100:100:1300]; clabel(CC,hh,vv,'FontSize',9,'Color',0.3*[1 1 1],'LabelSpacing',600)
%     contour(xx,yy,ps.z*reshape(s-b,gg.nI,gg.nJ),1*[1 1],'color','k','linewidth',1); % mark margin

    % add boundary
    plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5);
    
    % color bar
    tmp = get(gca,'position'); 
    cbar_size = [tmp(1)+tmp(3)-0.02 tmp(2)+tmp(4)*1/8 0.01 tmp(4)*6/8];
    cx2 = colorbar('peer',gca,'vertical','position',cbar_size,'YAxislocation','right');
    text(5,22,'H [ m ]','VerticalAlignment','bottom','HorizontalAlignment','right');

    cmap2=colormap(gca, parula);
    text(-0.11,0.5,'y [ km ]','HorizontalAlignment','center','units','normalized','rotation',90);
    set(gca,'TickDir','out','LineWidth',1.1); box on;  
    
axes(axe3)
    % velocity
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny;  
    zz = (ps.u_b*pd.ty)*reshape(aa.Ub,gg.nI,gg.nJ); xx(zz==0) = NaN;
    cax = [10^2 14000]; 
    
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis     
    hand = pcolor(xx,yy,zz); set(hand,'linestyle','none'); % shading interp;
    caxis(cax); 
    axis image;
    axis(axx);
    
    % add boundary
    plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5);
    
    % color bar
    tmp = get(gca,'position'); 
    cbar_size = [tmp(1)+tmp(3)-0.02 tmp(2)+tmp(4)*1/8 0.01 tmp(4)*6/8];
    cx3 = colorbar('peer',gca,'vertical','position',cbar_size,'YAxislocation','right');
    text(8,22,'U_{b} [ m y^{-1} ]','VerticalAlignment','bottom','HorizontalAlignment','right');

    load cmapu1; cmap3=cmap; colormap(gca, cmap3);
    if reversey, set(gca,'YDir','reverse'); end
    text(-0.11,0.5,'y [ km ]','HorizontalAlignment','center','units','normalized','rotation',90);
    text(0.5,-0.13/aspect,'x [ km ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized');
    set(gca,'TickDir','out','LineWidth',1.1); box on;

% surface input first calculate E (it's not automatically saved)
% for i_s=1:5   
%     load([fn,'/',int2four(i_s)]); 
%     %% extract new variables
%     if isfield(vv,'nbdy'), gg = nevis_label(gg,vv.nbdy); aa.phi = aa.phi_a(gg.nbdy); end % remask
%     aa = nevis_inputs(vv.t,aa,pp,gg,oo);
%     oo.evaluate_variables = 1; oo.evaluate_residual = 0; oo.evaluate_jacobian = 0; [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
%     vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
%     nevis_unpack(aa,vv2);
%     % save plotting variable E
%     E_i_s(:,i_s) = E; 
%     clear tt vv vv2
% end
% % average plotting variable over 6
% Emean = nanmean(E_i_s,2);
 
axes(axe7);
    % surface input
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny;  
    zz = (ps.m*10^3*pd.td)*reshape(E,gg.nI,gg.nJ); 
    cax = [0 150]; 
    
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis     
    hand = pcolor(xx,yy,zz); set(hand,'linestyle','none'); % shading interp;
    caxis(cax); 
    axis image;
    axis(axx);
    
    % color bar
    tmp = get(gca,'position'); 
    cbar_size = [tmp(1)+tmp(3)-0.02 tmp(2)+tmp(4)*1/8 0.01 tmp(4)*6/8];
    cx4 = colorbar('peer',gca,'vertical','position',cbar_size,'YAxislocation','right');
    text(8.5,22.5,'E [ mm d^{-1} ]','VerticalAlignment','bottom','HorizontalAlignment','right');
    time =  [tt.t]; DOY = time(end)*10;
    title(sprintf('DOY 2007 %6.2f',DOY))
    
%     % add moulins
%     if ~isfield(pp,'ni_m'), pp.ni_m = []; end
%     mscale = 100;  %amount to scale moulin input by to make size of moulin dot
%     x = (ps.x/10^3)*nx(pp.ni_m);
%     y = (ps.x/10^3)*ny(pp.ni_m);
%     hold on;
%     for i_m = 1:length(pp.ni_m),
%         if E(pp.ni_m(i_m))>0,
%             plot(x(i_m),y(i_m),'ko','Markersize',4+E(pp.ni_m(i_m))/mscale,'MarkerFaceColor',1*[1 1 1]); % mark moulins   
%         else
%             plot(x(i_m),y(i_m),'ko','Markersize',4,'MarkerFaceColor',0.8*[1 1 1]); % mark moulins  
%         end
%     end
            
    % add boundary
    plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5); 
    
    % add outlets
    plot((ps.x/10^3)*nx(gg.nbdy),(ps.x/10^3)*ny(gg.nbdy),'ko','markersize',4,'markerfacecolor','y');
    
    load cmapq2; cmap1=cmap; colormap(gca, cmap1);
    set(gca,'TickDir','out','LineWidth',1.1); box on;

axes(axe4);
    % basal melt input
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny;  
    zz = (ps.m*10^3*pd.td)*reshape(aa.m,gg.nI,gg.nJ); 
    cax = [0 40]; 
    
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis     
    hand = pcolor(xx,yy,zz); set(hand,'linestyle','none'); % shading interp;
    caxis(cax); 
    axis image;
    axis(axx);
    
    % color bar
    tmp = get(gca,'position'); 
    cbar_size = [tmp(1)+tmp(3)-0.02 tmp(2)+tmp(4)*1/8 0.01 tmp(4)*6/8];
    cx4 = colorbar('peer',gca,'vertical','position',cbar_size,'YAxislocation','right');
    text(9,22,'m [ mm d^{-1} ]','VerticalAlignment','bottom','HorizontalAlignment','right');

    % add moulins
    if ~isfield(pp,'ni_m'), pp.ni_m = []; end
    mscale = 100;  %amount to scale moulin input by to make size of moulin dot
    x = (ps.x/10^3)*nx(pp.ni_m);
    y = (ps.x/10^3)*ny(pp.ni_m);
    hold on;
    for i_m = 1:length(pp.ni_m),
        if E(pp.ni_m(i_m))>0,
            plot(x(i_m),y(i_m),'ko','Markersize',4+E(pp.ni_m(i_m))/mscale,'MarkerFaceColor',1*[1 1 1]); % mark moulins   
        else
            plot(x(i_m),y(i_m),'ko','Markersize',4,'MarkerFaceColor',0.8*[1 1 1]); % mark moulins  
        end
    end
    
    % add boundary
    plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5); 
    
    % add outlets
    plot((ps.x/10^3)*nx(gg.nbdy),(ps.x/10^3)*ny(gg.nbdy),'ko','markersize',4,'markerfacecolor','y');
    
    cmap4=load('cmapq2.mat');  colormap(gca, cmap4.cmap);
    set(gca,'TickDir','out','LineWidth',1.1); box on;
    
axes(axe5);
% discharge
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny;
    zz = ps.qs*reshape(qs+qe+qQ,gg.nI,gg.nJ);   
    cax = [10^(-5) 10^(0)]; 
    
    % logarithmic pcolor
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis
    hand = pcolor(xx,yy,log10(zz)); set(hand,'linestyle','none'); % shading interp;
    caxis(log10(cax)); 
    axis image;
    axis(axx);
    
    % color bar
    tmp = get(gca,'position'); 
    cbar_size = [tmp(1)+tmp(3)-0.02 tmp(2)+tmp(4)*1/8 0.01 tmp(4)*6/8];
    cx5 = colorbar('peer',gca,'vertical','position',cbar_size,'YAxislocation','right');
    text(8,22,'q [ m^2 s^{-1} ]','VerticalAlignment','bottom','HorizontalAlignment','right');
    tmp = get(cx5,'XTick');  labels = 10.^tmp;  set(cx5,'XTickLabel',labels);

    % add pressure contours
    hold on; hand = contour((ps.x/10^3)*nx,(ps.x/10^3)*ny,(ps.phi/10^6)*reshape(phi,gg.nI,gg.nJ),p_conts,'linecolor','r','linewidth',0.5);
    
%     % add moulins
%     if ~isfield(pp,'ni_m'), pp.ni_m = []; end
%     mscale = 100;  %amount to scale moulin input by to make size of moulin dot
%     x = (ps.x/10^3)*nx(pp.ni_m);
%     y = (ps.x/10^3)*ny(pp.ni_m);
%     hold on;
%     for i_m = 1:length(pp.ni_m),
%         if E(pp.ni_m(i_m))>0,
%             plot(x(i_m),y(i_m),'ko','Markersize',4+E(pp.ni_m(i_m))/mscale,'MarkerFaceColor',1*[1 1 1]); % mark moulins   
%         else
%             plot(x(i_m),y(i_m),'ko','Markersize',4,'MarkerFaceColor',0.8*[1 1 1]); % mark moulins  
%         end
%     end
    % add boundary
    plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5); 
    % add outlets
    plot((ps.x/10^3)*nx(gg.nbdy),(ps.x/10^3)*ny(gg.nbdy),'ko','markersize',4,'markerfacecolor','y');
    
    load cmapq2; cmap1=cmap; colormap(gca, cmap1);
    set(gca,'TickDir','out','LineWidth',1.1); box on;

axes(axe6)    
% effective pressure
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny;  
    zz = (ps.phi/10^6)*reshape(phi_0-phi,gg.nI,gg.nJ); 
    cax = [-1 1]; 
    
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis     
    hand = pcolor(xx,yy,zz); set(hand,'linestyle','none'); % shading interp;
    caxis(cax); 
    axis image;
    axis(axx);
    
    % color bar
    load cmapbluered; colormap(gca, cmap(end:-1:1,:));
    %cmap2=load('nevis/cmapbluered.mat'); cmap4=cmap2.cmap; colormap(gca, cmap4(end:-1:1,:));
    text(0.5,-0.13/aspect,'x [ km ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized');
   
    tmp = get(gca,'position'); 
    cbar_size = [tmp(1)+tmp(3)-0.02 tmp(2)+tmp(4)*1/8 0.01 tmp(4)*6/8];
    cx6 = colorbar('peer',gca,'vertical','position',cbar_size,'YAxislocation','right');
    text(7,22,'N [ MPa ]','VerticalAlignment','bottom','HorizontalAlignment','right');
    
%     % add moulins
%     if ~isfield(pp,'ni_m'), pp.ni_m = []; end
%     mscale = 100;  %amount to scale moulin input by to make size of moulin dot
%     x = (ps.x/10^3)*nx(pp.ni_m);
%     y = (ps.x/10^3)*ny(pp.ni_m);
%     hold on;
%     for i_m = 1:length(pp.ni_m),
%         if E(pp.ni_m(i_m))>0,
%             plot(x(i_m),y(i_m),'ko','Markersize',4+E(pp.ni_m(i_m))/mscale,'MarkerFaceColor',1*[1 1 1]); % mark moulins   
%         else
%             plot(x(i_m),y(i_m),'ko','Markersize',4,'MarkerFaceColor',0.8*[1 1 1]); % mark moulins  
%         end
%     end

    % add boundary
    plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5); 
    
    % add outlets
    plot((ps.x/10^3)*nx(gg.nbdy),(ps.x/10^3)*ny(gg.nbdy),'ko','markersize',4,'markerfacecolor','y');
    set(gca,'TickDir','out','LineWidth',1.1); box on;

axes(axe8) 
% sheet thickness
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny;  
    zz = (ps.h*100)*reshape(hs,gg.nI,gg.nJ); 
    cax = [0 50]; 
    
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis     
    hand = pcolor(xx,yy,zz); set(hand,'linestyle','none'); % shading interp;
    caxis(cax); 
    axis image;
    axis(axx);
    
    % color bar
    tmp = get(gca,'position'); 
    cbar_size = [tmp(1)+tmp(3)-0.02 tmp(2)+tmp(4)*1/8 0.01 tmp(4)*6/8];
    cx8 = colorbar('peer',gca,'vertical','position',cbar_size,'YAxislocation','right');
    text(7.75,22,'h_s [ cm ]','VerticalAlignment','bottom','HorizontalAlignment','right');
    
    load cmapq2; cmap3=cmap; colormap(gca, cmap3);
    if reversey, set(gca,'YDir','reverse'); end

%     % add moulins
%     if ~isfield(pp,'ni_m'), pp.ni_m = []; end
%     mscale = 100;  %amount to scale moulin input by to make size of moulin dot
%     x = (ps.x/10^3)*nx(pp.ni_m);
%     y = (ps.x/10^3)*ny(pp.ni_m);
%     hold on;
%     for i_m = 1:length(pp.ni_m),
%         if E(pp.ni_m(i_m))>0,
%             plot(x(i_m),y(i_m),'ko','Markersize',4+E(pp.ni_m(i_m))/mscale,'MarkerFaceColor',1*[1 1 1]); % mark moulins   
%         else
%             plot(x(i_m),y(i_m),'ko','Markersize',4,'MarkerFaceColor',0.8*[1 1 1]); % mark moulins  
%         end
%     end
    
    % add boundary
    plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5); 
    % add outlets
    plot((ps.x/10^3)*nx(gg.nbdy),(ps.x/10^3)*ny(gg.nbdy),'ko','markersize',4,'markerfacecolor','y');
    set(gca,'TickDir','out','LineWidth',1.1); box on;
    
axes(axe9);
 % channel volume, converted to sheet thickness
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny;  
    zz = (ps.h/10)*reshape(VS./(gg.Dx.*gg.Dy),gg.nI,gg.nJ); 
    cax = [0 0.001]; 
    
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis     
    hand = pcolor(xx,yy,zz); set(hand,'linestyle','none'); % shading interp;
    caxis(cax); 
    axis image;
    axis(axx);
    
    % color bar
    tmp = get(gca,'position'); 
    cbar_size = [tmp(1)+tmp(3)-0.02 tmp(2)+tmp(4)*1/8 0.01 tmp(4)*6/8];
    cx9 = colorbar('peer',gca,'vertical','position',cbar_size,'YAxislocation','right');
    text(8,22.75,'h_{VS} [ cm ]','VerticalAlignment','bottom','HorizontalAlignment','right');
    
    load cmapq2; cmap2=cmap; colormap(gca, cmap2);
    text(0.5,-0.12/aspect,'x [ km ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized');
    
    % add pressure contours
    hold on; hand = contour((ps.x/10^3)*nx,(ps.x/10^3)*ny,(ps.phi/10^6)*reshape(phi,gg.nI,gg.nJ),p_conts,'linecolor','r','linewidth',0.5);
       
%     % add moulins
%     if ~isfield(pp,'ni_m'), pp.ni_m = []; end
%     mscale = 100;  %amount to scale moulin input by to make size of moulin dot
%     x = (ps.x/10^3)*nx(pp.ni_m);
%     y = (ps.x/10^3)*ny(pp.ni_m);
%     hold on;
%     for i_m = 1:length(pp.ni_m),
%         if E(pp.ni_m(i_m))>0,
%             plot(x(i_m),y(i_m),'ko','Markersize',4+E(pp.ni_m(i_m))/mscale,'MarkerFaceColor',1*[1 1 1]); % mark moulins   
%         else
%             plot(x(i_m),y(i_m),'ko','Markersize',4,'MarkerFaceColor',0.8*[1 1 1]); % mark moulins  
%         end
%     end

    % add boundary
    plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5);  
    
    % add outlets
    plot((ps.x/10^3)*nx(gg.nbdy),(ps.x/10^3)*ny(gg.nbdy),'ko','markersize',4,'markerfacecolor','y');   
    set(gca,'TickDir','out','LineWidth',1.1); box on;
    
print(gcf,'-dpng','-r500',sprintf(['panel9_h22222_ubspatial_restart_lakeramp2_%.png',DOY]));