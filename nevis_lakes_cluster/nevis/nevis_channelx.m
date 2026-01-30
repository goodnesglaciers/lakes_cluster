%function nevis_channelx(vv,aa,pp,ps,gg,oo) % [uncomment to make a function ]
% nevis_channelx(vv,aa,pp,ps,gg,oo)
% comment first line to use as script taking inputs from current workspace
% plots direrctions of channels (Sx, Sy, Sr, Sc) in a four-panel figure
% with channel cross-sectional area
% 
% 17 July 2016: taken from nevis_plot for regional channel cross sections (LAS) 

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
z_conts = -1000:50:5000;
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
if ~isempty(gg.n1),
x_out = gg.nx(gg.n1); y_out = gg.ny(gg.n1);
tmp = nevis_orderboundary(x_out,y_out); x_out = x_out(tmp); y_out = y_out(tmp); % reorder to follow boundary
else x_out = []; y_out = [];
end

%% axes limits
axx = (ps.x/10^3)*[gg.xl gg.xr gg.yb gg.yt];

%% set up figure
figure; clf; set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',[1 1 20 20]);
set(0,'DefaultAxesFontSize',10);
set(0,'DefaultTextFontSize',10);
axes_size = [0.1 0.1 0.4 0.4];

axe1 = axes('Position',[0.08 0.65 0.4 0.3],'Box','on','NextPlot','add'); 
axe2 = axes('Position',[0.55 0.65 0.4 0.3],'Box','on','NextPlot','add'); 
axe3 = axes('Position',[0.08 0.33 0.4 0.3],'Box','on','NextPlot','add'); 
axe4 = axes('Position',[0.55 0.33 0.4 0.3],'Box','on','NextPlot','add'); 

axe5 = axes('Position',[0.08 0.02 0.4 0.3],'Box','on','NextPlot','add');

%adjust axes positions to account for aspect ratio of figure
aspect = (axx(4)-axx(3))/(axx(2)-axx(1)); % aspect of axes limits
pos = axes_size; tmp = get(gcf,'position'); sca = tmp(4)/tmp(3);
if aspect>(pos(4)/pos(3))*sca, 
    pos = pos + [pos(4)*sca*(1-1/aspect)/2 0 -pos(3)+pos(4)*sca*(1/aspect) 0]; % constrained by y
else
    pos = pos + [0 pos(3)/sca*(1-aspect)/2 0 -pos(4)+pos(3)/sca*aspect]; % constrained by x
end 
axes_size = pos;

%% plot

% NOTES LAS
% % x direction
% imagesc(ps.x*gg.ex(:,1),ps.x*gg.ey(1,:),reshape(log10(ps.S*vv.Sx),gg.eI,gg.eJ)'
% % other direction diagonoal
% imagesc(ps.x*gg.cx(:,1),ps.x*gg.cy(1,:),reshape(ps.S*vv.Sr,gg.cI,gg.cJ)');

% FOR AREA PLOTS
%     % channel volume, converted to sheet thickness
%     xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny;  
%     zz = (ps.h/10)*reshape(VS./(gg.Dx.*gg.Dy),gg.nI,gg.nJ); 
%     cax = [0 0.001]; 

if channelx
axes(axe1);
    % channelx x-direction (west)
    xx = (ps.x/10^3)*ex; yy = (ps.x/10^3)*ey;
    zz = reshape((ps.S*vv.Sx)*10,gg.eI,gg.eJ);   
    cax = [10^(-2) 10^(2)]; 
        
    % logarithmic pcolor
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis
    hand = pcolor(xx,yy,log10(zz)); set(hand,'linestyle','none'); % shading interp;
    caxis(log10(cax)); 
    axis image;
    axis(axx);
    
    % color bar
    tmp = get(gca,'position'); cbar_size = [tmp(3)*3/16 0.01  tmp(3)*5/8 0.02];
    cx = colorbar('peer',gca,'horizontal','position',[0.35 0.95 0.4 0.02],'XAxislocation','top');
    tmp = get(cx,'XTick'); labels = 10.^tmp; set(cx,'XTickLabel',labels);
    text(0.55,1.2,'h_{S} [ cm ]','units','normalized','VerticalAlignment','bottom','HorizontalAlignment','right');
%    cx.Label.String = 'q [ m^2 s^{-1} ]'; cx.Label.Units = 'normalized'; cx.Label.Position = [-0.05 0]; cx.Label.VerticalAlignment = 'bottom'; cx.Label.HorizontalAlignment = 'right';
    
    % add pressure contours
    hold on; hand = contour((ps.x/10^3)*nx,(ps.x/10^3)*ny,(ps.phi/10^6)*reshape(phi,gg.nI,gg.nJ),p_conts,'linecolor','r','linewidth',0.5);
    
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
    
    % add outlets
    plot((ps.x/10^3)*nx(gg.nbdy),(ps.x/10^3)*ny(gg.nbdy),'ko','markersize',4,'markerfacecolor','y');
    
    % add boundary
    plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5); 
    
    load cmapq2; colormap(cmap);
    title('Sx (E-W)')
    if reversey, set(gca,'YDir','reverse'); end
    text(-0.12,0.5,'y [ km ]','HorizontalAlignment','center','units','normalized','rotation',90);
    text(0.5,-0.12/aspect,'x [ km ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized');
    if save_plot, print(gcf,'-depsc2',[fn,'_discharge']); end

axes(axe2);
    % channelx y-direction (north)
    xx = (ps.x/10^3)*fx; yy = (ps.x/10^3)*fy;
    zz = reshape((ps.S*vv.Sy)*10,gg.fI,gg.fJ);   
    cax = [10^(-2) 10^(2)];  
    
    % logarithmic pcolor
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis
    hand = pcolor(xx,yy,log10(zz)); set(hand,'linestyle','none'); % shading interp;
    caxis(log10(cax)); 
    axis image;
    axis(axx);
    
%     % color bar
%     tmp = get(gca,'position'); cbar_size = [tmp(3)*3/16 0.01  tmp(3)*5/8 0.02];
%     cx = colorbar('peer',gca,'horizontal','position',[tmp(1) tmp(2)+tmp(4) 0 0]+cbar_size,'XAxislocation','top');
%     tmp = get(cx,'XTick'); labels = 10.^tmp; set(cx,'XTickLabel',labels);
%     text(0.1,1.05,'q [ m^2 s^{-1} ]','units','normalized','VerticalAlignment','bottom','HorizontalAlignment','right');
%     cx.Label.String = 'q [ m^2 s^{-1} ]'; cx.Label.Units = 'normalized'; cx.Label.Position = [-0.05 0]; cx.Label.VerticalAlignment = 'bottom'; cx.Label.HorizontalAlignment = 'right';
    
    % add pressure contours
    hold on; hand = contour((ps.x/10^3)*nx,(ps.x/10^3)*ny,(ps.phi/10^6)*reshape(phi,gg.nI,gg.nJ),p_conts,'linecolor','r','linewidth',0.5);
    
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
    
    % add outlets
    plot((ps.x/10^3)*nx(gg.nbdy),(ps.x/10^3)*ny(gg.nbdy),'ko','markersize',4,'markerfacecolor','y');
    
    % add boundary
    plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5); 
    
    load cmapq2; colormap(cmap);
    title('Sy (N-S)')
    if reversey, set(gca,'YDir','reverse'); end
    text(-0.12,0.5,'y [ km ]','HorizontalAlignment','center','units','normalized','rotation',90);
    text(0.5,-0.12/aspect,'x [ km ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized');
    if save_plot, print(gcf,'-depsc2',[fn,'_discharge']); end

axes(axe3)    
    % channelx Sr diagonal (north west)
    xx = (ps.x/10^3)*gg.cx; yy = (ps.x/10^3)*gg.cy;
    zz = reshape((ps.S*vv.Sr)*10,gg.cI,gg.cJ);    
    cax = [10^(-2) 10^(2)];  
    
    % logarithmic pcolor
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis
    hand = pcolor(xx,yy,log10(zz)); set(hand,'linestyle','none'); % shading interp;
    caxis(log10(cax)); 
    axis image;
    axis(axx);
    
%     % color bar
%     tmp = get(gca,'position'); cbar_size = [tmp(3)*3/16 0.01  tmp(3)*5/8 0.02];
%     cx = colorbar('peer',gca,'horizontal','position',[tmp(1) tmp(2)+tmp(4) 0 0]+cbar_size,'XAxislocation','top');
%     tmp = get(cx,'XTick'); labels = 10.^tmp; set(cx,'XTickLabel',labels);
%     text(0.1,1.05,'q [ m^2 s^{-1} ]','units','normalized','VerticalAlignment','bottom','HorizontalAlignment','right');
%     cx.Label.String = 'q [ m^2 s^{-1} ]'; cx.Label.Units = 'normalized'; cx.Label.Position = [-0.05 0]; cx.Label.VerticalAlignment = 'bottom'; cx.Label.HorizontalAlignment = 'right';
    
    % add pressure contours
    hold on; hand = contour((ps.x/10^3)*nx,(ps.x/10^3)*ny,(ps.phi/10^6)*reshape(phi,gg.nI,gg.nJ),p_conts,'linecolor','r','linewidth',0.5);
    
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
    
    % add outlets
    plot((ps.x/10^3)*nx(gg.nbdy),(ps.x/10^3)*ny(gg.nbdy),'ko','markersize',4,'markerfacecolor','y');
    
    % add boundary
    plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5); 
    
    load cmapq2; colormap(cmap);
    title('Sr (NW-SE)')
    if reversey, set(gca,'YDir','reverse'); end
    text(-0.12,0.5,'y [ km ]','HorizontalAlignment','center','units','normalized','rotation',90);
    text(0.5,-0.12/aspect,'x [ km ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized');
    if save_plot, print(gcf,'-depsc2',[fn,'_discharge']); end
    
axes(axe4)    
    % channelx Ss diagonal (south west)
    xx = (ps.x/10^3)*gg.cx; yy = (ps.x/10^3)*gg.cy;
    zz = reshape((ps.S*vv.Ss)*10,gg.cI,gg.cJ); 
    cax = [10^(-2) 10^(2)];  
    
    % logarithmic pcolor
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis
    hand = pcolor(xx,yy,log10(zz)); set(hand,'linestyle','none'); % shading interp;
    caxis(log10(cax)); 
    axis image;
    axis(axx);
    
%     % color bar
%     tmp = get(gca,'position'); cbar_size = [tmp(3)*3/16 0.01  tmp(3)*5/8 0.02];
%     cx = colorbar('peer',gca,'horizontal','position',[tmp(1) tmp(2)+tmp(4) 0 0]+cbar_size,'XAxislocation','top');
%     tmp = get(cx,'XTick'); labels = 10.^tmp; set(cx,'XTickLabel',labels);
%     text(0.1,1.05,'q [ m^2 s^{-1} ]','units','normalized','VerticalAlignment','bottom','HorizontalAlignment','right');
%     cx.Label.String = 'q [ m^2 s^{-1} ]'; cx.Label.Units = 'normalized'; cx.Label.Position = [-0.05 0]; cx.Label.VerticalAlignment = 'bottom'; cx.Label.HorizontalAlignment = 'right';
    
    % add pressure contours
    hold on; hand = contour((ps.x/10^3)*nx,(ps.x/10^3)*ny,(ps.phi/10^6)*reshape(phi,gg.nI,gg.nJ),p_conts,'linecolor','r','linewidth',0.5);
  
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
    
    % add outlets
    plot((ps.x/10^3)*nx(gg.nbdy),(ps.x/10^3)*ny(gg.nbdy),'ko','markersize',4,'markerfacecolor','y');
    
    % add boundary
    plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5); 
    
    load cmapq2; colormap(cmap);
    title('Sr (NE-SW)')
    if reversey, set(gca,'YDir','reverse'); end
    text(-0.12,0.5,'y [ km ]','HorizontalAlignment','center','units','normalized','rotation',90);
    text(0.5,-0.12/aspect,'x [ km ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized');
    if save_plot, print(gcf,'-depsc2',[fn,'_discharge']); end
end

if area
% clf; axes('position',axes_size);
axes(axe5) 
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
    tmp = get(gca,'position'); cbar_size = [tmp(3)*3/16 0.01  tmp(3)*5/8 0.02];
    cx = colorbar('peer',gca,'horizontal','position',[tmp(1)+0.4 tmp(2)+tmp(4)-0.2 0 0]+cbar_size,'xAxislocation','top');
    text(1.3,0.5,'h_S [ cm ]','units','normalized');

    load cmapq2; colormap(cmap);
    if reversey, set(gca,'YDir','reverse'); end
    text(-0.12,0.5,'y [ km ]','HorizontalAlignment','center','units','normalized','rotation',90);
    text(0.5,-0.12/aspect,'x [ km ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized');
    if save_plot, print(gcf,'-depsc2',[fn,'_area']); end
end


