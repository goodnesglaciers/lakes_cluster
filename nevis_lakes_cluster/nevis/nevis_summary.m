% function nevis_plot_summary(tt,vv,aa,pp,ps,gg,oo) % [uncomment to make a function ]
% plot summary of time series in struct tt
%
% 21 August 2014

    t = (ps.t/(24*60*60))*[tt.t];
    Q_in = ps.Q*[tt.Q_in];
    Q_out = ps.Q*[tt.Q_out];
    m = (ps.m*ps.x^2)*[tt.m];
    E = (ps.m*ps.x^2)*[tt.E];

    phi = (ps.phi/10^6)*[tt.phi];
    N = (ps.phi/10^6)*[tt.N];
    hs = ps.x^2*ps.h*[tt.hs];
    he = ps.x^2*ps.h*[tt.he];
    S = ps.x*ps.S*[tt.S];
    A = ps.x^2*sum(gg.Dx.*gg.Dy);
    if isfield(tt,'pts_phi')    
    pts_phi = (ps.phi/10^6)*[tt.pts_phi];
    pts_hs = ps.hs*[tt.pts_hs];
    pts_N = (ps.phi/10^6)*(aa.phi_0(oo.pts_ni)*[tt.t].^0 - [tt.pts_phi]);
    pts_pw = (ps.phi/10^6)*([tt.pts_phi]-aa.phi_a(oo.pts_ni)*[tt.t].^0);
    pts_prat = ([tt.pts_phi]-aa.phi_a(oo.pts_ni)*[tt.t].^0)./(aa.phi_0(oo.pts_ni)*[tt.t].^0-aa.phi_a(oo.pts_ni)*[tt.t].^0);
    end

    figure(2); 
    ax(1) = subplot(3,1,1);
        plot(t,Q_in,'k',t,Q_out,'r',t,m,'b--',t,E,'b');
        xlabel('t [ d ]');
        ylabel('Q [ m^3/s ]');
        legend('Q_in','Q_out','m','E');
        ylim([0 max(E)])

    ax(2) = subplot(3,1,2); 
        if isfield(tt,'pts_phi')
        plot(t,pts_N,'-');
        end
        hold on;
        plot(t,N,'-','color',0.6*[1 1 1]);
        xlabel('t [ d ]');
        ylabel('N [ MPa ]');
        
    ax(3) = subplot(3,1,3); 
        plot(t,hs./A,'b-',t,he./A,'k-',t,1e2*S./A,'r-');
        xlabel('t [ d ]');
        ylabel('h [ m ]');
        legend('hcav (hs)','hel','S')