%function nevis_summary2(tt,vv,aa,pp,ps,pd,gg,oo) % [uncomment to make a function ]
% plot summary of time series in struct tt
%
% 31 October 2015 (Laura)

    cd(oo.fn)
    fnames = dir('*.mat');
    load(fnames(end).name)
    % figure name
    str = sprintf('Run %s',[oo.root,oo.fn]);
    str_save = sprintf('%s_summary2',[oo.root,oo.fn]);
    cd ../
    
    oo.evaluate_variables = 1; 
    [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
    vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
  
    % scaling back to dimensional quantities
    t = (ps.t/(24*60*60))*[tt.t];
    Q_in = ps.Q*[tt.Q_in];
    Q_out = ps.Q*[tt.Q_out];
    m = (ps.m*ps.x^2)*[tt.m];
    E = (ps.m*ps.x^2)*[tt.E];

    phi = (ps.phi/10^6)*[tt.phi];
    N = (ps.phi/10^6)*[tt.N];
    pw = (ps.phi/10^6)*([tt.phi]-aa.phi_a(1));
    hs = ps.hs*[tt.hs];
    he = ps.he*[tt.he];
    S = ps.x*ps.S*[tt.S];
    A = ps.x^2*sum(gg.Dx.*gg.Dy);
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

    fig2=figure(2); 
    %Q
        ax(1) = subplot(5,1,1);
        plot(t,Q_in,'k',t,Q_out,'r',t,m,'b--',t,E,'b');
        %xlabel('t [ d ]');
        ylabel('Q [ m^3/s ]');
        legend('Q_in','Q_out','m','E');
        title(str);
   %phi
        ax(2) = subplot(5,1,2); 
        if isfield(tt,'pts_phi')
        plot(t,pts_N,'-');
        end
        hold on;
        plot(t,N,'-','color',0.6*[1 1 1]);
        %xlabel('t [ d ]');
        ylabel('N [ MPa ]');
        %legend('Lake','[-4 0]','[-8 0]','avg');      
        
   %S
        ax(3) = subplot(5,1,3); 
        if isfield(tt,'pts_phi')
        plot(t,pts_S./pts_A,'-');
        end
        hold on;
        %plot(t,S./A,'-','color',0.6*[1 1 1]);
        %%plot(t,hs./A,'b-',t,he./A,'k-',t,1e2*S./A,'r-');
        %xlabel('t [ d ]');
        ylabel('S [ m^{2} ?]');
        %legend('[0 0]','[-4 0]','[-8 0]','avg'); 
        
    %hcav
        ax(4) = subplot(5,1,4); 
        if isfield(tt,'pts_phi')
        plot(t,pts_hs,'-');
        end
        hold on;
        %plot(t,hs./A,'-','color',0.6*[1 1 1]);   
        %xlabel('t [ d ]');
        ylabel('h_{cav} [ m ]');
        %legend('[0 0]','[-4 0]','[-8 0]','avg')
   
   %he
   ax(5) = subplot(5,1,5);
        if isfield(tt,'pts_phi')
        plot(t,pts_he,'-');
        end
        hold on;
        plot(t,he,'-','color',0.6*[1 1 1]);
        xlabel('t [ d ]');
        ylabel('h_{el} [ m ]');
        legend('[0 0]','[-4 0]','[-8 0]','avg')

%    %Pw 
%    ax(5) = subplot(5,1,5);
%         if isfield(tt,'pts_phi')
%         plot(t,pts_pw,'-');
%         end
%         hold on;
%         plot(t,pw,'-','color',0.6*[1 1 1]);
%         xlabel('t [ d ]');
%         ylabel('P_{w} [ MPa ]');
%         legend('[0 0]','[-4 0]','[-8 0]','avg')

set(fig2,'PaperPositionMode', 'auto')
print(fig2,str_save,'-depsc2')

%% nevis_summary3 transect plots

x = [-1 1];             % start end X
y = [0 0];                  % start end Y
n = 100;                    % number of points in transect
% N (from nevis_animate)
%z_N = (ps.phi/10^6)*reshape(aa.phi_0-vv.phi,gg.nI,gg.nJ); % variable of interest for transect at nodes
z_N = (ps.phi/10^6)*reshape(vv2.N,gg.nI,gg.nJ);
[x_tran,y_tran,s_tran,N_tran] = nevis_phi_transect(x,y,n,z_N,gg);

% Q discharge (from nevis_animate)
z_Q = ps.qs*reshape(vv2.qs+vv2.qe+vv2.qQ,gg.nI,gg.nJ);  
[x_tran,y_tran,s_tran,Q_tran] = nevis_phi_transect(x,y,n,z_Q,gg);

% Q mean channel discharge on nodes, scaled with ps.Q
z_Q2 = ps.Q*reshape(vv2.Q,gg.nI,gg.nJ);  
[x_tran,y_tran,s_tran,Q2_tran] = nevis_phi_transect(x,y,n,z_Q2,gg);

% S (volume of channels at nodes, scaled with ps.h*ps.x^2)
z_S = (ps.h*ps.x^2)*reshape(vv2.VS,gg.nI,gg.nJ);
[x_tran,y_tran,s_tran,S_tran] = nevis_phi_transect(x,y,n,z_S,gg);

% he
z_he = (ps.h)*reshape(vv2.he,gg.nI,gg.nJ);
[x_tran,y_tran,s_tran,he_tran] = nevis_phi_transect(x,y,n,z_he,gg);
% hs
z_hs = (ps.h)*reshape(vv.hs,gg.nI,gg.nJ);
[x_tran,y_tran,s_tran,hs_tran] = nevis_phi_transect(x,y,n,z_hs,gg);

fig3=figure(3)

    ax(1) = subplot(2,1,1)
    [ax1,h1, h2] =  plotyy(ps.x*x_tran/1e3,N_tran,ps.x*x_tran/1e3,Q_tran,'line','line')
    hold all
    plot(ps.x*x_tran/1e3,Q2_tran,'Parent',ax1(2))
    ylabel(ax1(1),'N [ MPa ]') % left y-axis
    ylabel(ax1(2),'Q [ m^{3}/s ]') % right y-axis
    %xlabel(ax1(1),' x  (km)');
    set(h1,'LineWidth',1.1,'color','b')
    set(h2,'LineWidth',1.1,'LineStyle','-')
    
    ax(2) = subplot(2,1,2)
    [ax2,h1, h2] =  plotyy(ps.x*x_tran/1e3,hs_tran,...
        ps.x*x_tran/1e3,S_tran,'line','line')
    hold on
    %plot(ps.x*x_tran/1e3,he_tran,'Parent',ax2(2));
    ylabel(ax2(1),'h_{cav} and h_{el} [ m ]') % left y-axis
    ylabel(ax2(2),'Volume at nodes [ m^{3} ]') % right y-axis
    xlabel(ax2(1),' x  (km)');
    set(h1,'LineWidth',1.1,'color','b')
    set(h2,'LineWidth',1.1,'LineStyle','--','color','r')
    

str_save3 = sprintf('%s_summary3',[oo.root,oo.fn]);   
set(fig3,'PaperPositionMode', 'auto')
print(fig3,str_save3,'-depsc2')