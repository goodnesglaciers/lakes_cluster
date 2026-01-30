function nevis_tranect_plot(fn,tis,transectx,transecty,n_tran);      %tt,vv,aa,pp,ps,pd,gg,oo) % [uncomment to make a function ]
% plot summary of transect time series from tt for one run
% 2 November 2015 (Laura)
% Inputs:
%   fn            filename
%   tis           time steps for plotting
%   n_tran        number of points in transect

%% load initial timestep
if exist([fn,'/0000']), load([fn,'/0000']); else load(fn); end
load([fn,'/0001']);
str = sprintf('%s.mat',[fn]); load(str)
str_save = sprintf('%s_transect',[oo.root,oo.fn]);
if isfield(oo,'reversey'), reversey = oo.reversey; else reversey = 0; end

%% set up figure
fig2=figure(2); clf; set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',[5 2 20 20]);
set(0,'DefaultAxesFontSize',10);
set(0,'DefaultTextFontSize',10);
axes_size = [0.1 0.1 0.8 0.8];
axe1 = axes('Position',[ 0.1 0.80 0.8 0.16 ],'Box','on','NextPlot','add'); 
axe2 = axes('Position',[ 0.1 0.58 0.8 0.16 ],'Box','on','NextPlot','add'); 
axe3 = axes('Position',[ 0.1 0.405 0.8 0.16 ],'Box','on','NextPlot','add');
axe4 = axes('Position',[ 0.1 0.235 0.8 0.16 ],'Box','on','NextPlot','add');
axe5 = axes('Position',[ 0.1 0.06 0.8 0.16 ],'Box','on','NextPlot','add');

%% plot timesteps
for i_t = 1:length(tis)
    
%% load timestep
load([fn,'/',int2four(tis(i_t))]);

%% extract initial variables
if isfield(vv,'nbdy'), gg = nevis_label(gg,vv.nbdy); aa.phi = aa.phi_a(gg.nbdy); end % re-mask
oo.evaluate_variables = 1; oo.evaluate_residual = 1; oo.evaluate_jacobian = 0; [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
%nevis_unpack(aa,gg);

%get rid of points outside domain
nx(gg.nout) = NaN;
ex(gg.eout) = NaN;
fx(gg.fout) = NaN;
cx(gg.cout) = NaN;
    
%boundary curve
x_out = gg.nx(gg.n1); y_out = gg.ny(gg.n1);
if ~isempty(x_out), tmp = nevis_orderboundary(x_out,y_out); x_out = x_out(tmp); y_out = y_out(tmp); end % reorder to follow boundary

% %% extract new variables
% if isfield(vv,'nbdy'), gg = nevis_label(gg,vv.nbdy); aa.phi = aa.phi_a(gg.nbdy); end % remask
% aa = nevis_inputs(vv.t,aa,pp,gg,oo);
% oo.evaluate_variables = 1; oo.evaluate_residual = 1; oo.evaluate_jacobian = 0; [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
% vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
% % nevis_unpack(aa,vv2);

    t = (ps.t/(24*60*60))*[tt.t];
    Q_in = ps.Q*[tt.Q_in];
    Q_out = ps.Q*[tt.Q_out];
    Q_outQ = ps.Q*[tt.Q_outQ];
    Q_outq = ps.Q*[tt.Q_outq];
    m = (ps.m*ps.x^2)*[tt.m];
    E = (ps.m*ps.x^2)*[tt.E];

    phi = (ps.phi/10^6)*[tt.phi];
    N = (ps.phi/10^6)*[tt.N];
    pw = (ps.phi/10^6)*([tt.phi]-aa.phi_a(1));
    hs = ps.x^2*ps.h*[tt.hs];               % total cavity sheet volume scaled
    he = ps.x^2*ps.h*[tt.he];               % total hel volume scaled
    S = ps.x*ps.S*[tt.S];                   % total channel vilume scaled
    A = ps.x^2*sum(gg.Dx.*gg.Dy);           % area of the domain
    total_water = S + hs + he;              % total volume of water in domain (hcav+hel+channel)
    
    if isfield(tt,'pts_phi')    
    pts_phi = (ps.phi/10^6)*[tt.pts_phi];
    pts_hs = ps.hs*[tt.pts_hs];
    pts_he = ps.he*[tt.pts_he];
    pts_N = (ps.phi/10^6)*(aa.phi_0(oo.pts_ni)*[tt.t].^0 - [tt.pts_phi]);
    pts_pw = (ps.phi/10^6)*([tt.pts_phi]-aa.phi_a(oo.pts_ni)*[tt.t].^0);
    pts_prat = ([tt.pts_phi]-aa.phi_a(oo.pts_ni)*[tt.t].^0)./(aa.phi_0(oo.pts_ni)*[tt.t].^0-aa.phi_a(oo.pts_ni)*[tt.t].^0);
    pts_S = ps.x*ps.S*([tt.pts_Sx]+[tt.pts_Sy]+[tt.pts_Sr]+[tt.pts_Ss]);
    pts_A = ps.x^2*4*(gg.Dx(1)*gg.Dy(1));
    end

% nevis_var_transext_xe (obtain variable on edges in x-direction)        
xt = [-1 1];        % start end X
yt = [0 0];         % start end Y
nt = n_tran;        % number of points in transect
% Q discharge (from nevis_animate)
%z_Q = ps.qs*reshape(vv2.qs+vv2.qe+vv2.qQ,gg.nI,gg.nJ);   % magnitude on nodes
%[x_tran,y_tran,s_tran,Q_tran] = nevis_var_transect_xe(xt,yt,nt,z_Q,gg);
z_Qx = ps.qs*reshape(vv2.Qx,gg.eI, gg.eJ);               % Q on edges e (x)
[x_tran,y_tran,s_tran,Qx_tran] = nevis_var_transect_xe(xt,yt,nt,z_Qx,gg);

% x and y GPS points
cc=jet(length(tis));

axes(axe2)
plot(ps.x*x_tran/1e3,Qx_tran,'-','Color',cc(i_t,:))
set(axe2,'xticklabel',[]); ylabel('Q_{x}  [ m^{3}/s ]');
% cb=colorbar('position',[0.95 0.75 0.2 0.03],'xlabel','i_t');
hold all
axes(axe3)
plot(transectx/1e3,pts_phi(:,tis(i_t)),'-','Color',cc(i_t,:))
set(axe3,'xticklabel',[]); ylabel('Phi [ MPa ]');
hold all
axes(axe4)
plot(transectx/1e3,pts_hs(:,tis(i_t)),'-','Color',cc(i_t,:))
set(axe4,'xticklabel',[]); ylabel('h_{cav} [ m ]');
hold all
axes(axe5)
plot(transectx/1e3,pts_he(:,tis(i_t)),'-','Color',cc(i_t,:))
hold on
str_var = sprintf('S_{0} = %.2f;  l_{r} = %.2f;  h_{r} = %.2f;  k_{c} = %.2f',pd.S_0,pd.l_r,pd.h_r,pd.k_c);
text(-9.5, 1, str_var);
xlabel(' Distance along transect [ km ]'); ylabel('h_{el} [ m ]')
hold all

end
axes(axe1)
plot(t, Q_out,'k','LineWidth',1.1); %total Q out
hold on
plot(t, Q_in,'Color',[0.7 0.7 0.7],'LineWidth',1.1); % Q_in
plot(t, Q_outQ,'r--','LineWidth',1.1); %channel
plot(t, Q_outq,'b--','LineWidth',1.1); %sheet
legend('Q_{out}','Q_{in}','Q_{channel}','Q_{sheet}','Location','NorthEast');
xlabel('t [ d ]'); ylabel('Q  [ m^{3}/s ]');
str_title = sprintf('%s transect timeseries',[oo.root,oo.fn]); title(str_title);

clear tt vv vv2

set(fig2,'PaperPositionMode', 'auto')
print(fig2,str_save,'-depsc2')

fig3=figure(3); % Volume conservation
clf
subplot(4,1,1)
plot(t, (hs)/1e9, 'k')
ylabel('Volume (hs) [ km^{3} ]'); xlabel('t [ d ]');

subplot(4,1,2)
plot(t, (he)/1e9, 'k')
ylabel('Volume (he) [ km^{3} ]'); xlabel('t [ d ]');

subplot(4,1,3)
plot(t, (S)/1e9, 'k')
ylabel('Volume (S) [ km^{3} ]'); xlabel('t [ d ]');

subplot(4,1,4)
plot(t, (S+hs+he)/1e9, 'k')
hold on
plot(t, (he+hs)/1e9,'b')
legend('S+hs+he','hs+he','Location','NorthWest');
ylabel('Volume (S+hs+he) [ km^{3} ]'); xlabel('t [ d ]');

figure(4)
plot(t, he/1e9, 'b')
hold on
plot(t, hs/1e9, 'k')
plot(t, S/1e9, 'r')
legend('V_{hel}','V_{hcav}','V_{channels}','Location','NorthWest');
ylabel('Volume [ km^{3} ]'); xlabel('t [ d ]');

str_save2 = sprintf('%s_transect2',[oo.root,oo.fn]);
set(fig3,'PaperPositionMode', 'auto')
print(fig3,str_save2,'-depsc2')
% 
% figure(5)
% plot(t, (he+hs)/1e9)
end
