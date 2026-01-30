%function nevis_summary2(tt,vv,aa,pp,ps,pd,gg,oo) % [uncomment to make a function ]
% plot summary of time series in struct tt
%
% 8 June 2017 (Laura Stevens)

%     cd(oo.fn)
%     fnames = dir('*.mat');
%     load(fnames(end).name)
%     % figure name
%     str = sprintf('Run %s',[oo.root,oo.fn]);
%     str_save = sprintf('%s_summary2',[oo.root,oo.fn]);
%     cd ../
    
    oo.evaluate_variables = 1; 
    [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
    vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
  
    % scaling back to dimensional quantities
    t = (ps.t/(24*60*60))*[tt.t];
    Q_in = ps.Q*[tt.Q_in];
    Q_out = ps.Q*[tt.Q_out];
    m = (ps.m*ps.x^2)*[tt.m];
    E = (ps.m*ps.x^2)*[tt.E];

    phi = (ps.phi/10^6)*[tt.phi]; %MPa
    N = (ps.phi/10^6)*[tt.N];      % MPa
    pw = (ps.phi/10^6)*([tt.phi]-aa.phi_a(1)); % MPa
    hs = ps.hs*(ps.x^2)*[tt.hs]; 
    he = ps.he*(ps.x^2)*[tt.he];
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
%%
    % save workspace to keep variables for plotting
    save varh22222_ubspatial_R67_lakerampM_4tiles_Ks100H.mat
    str='varh22222_ubspatial_R67_lakerampM_4tiles_Ks100H'
    save tt_ubspatial_R67_lakerampM_4tiles_Ks100H.mat tt
 
    
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
        legend('lake','29','28','27','26','39','38','36','35');
        
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
        %plot(t,he,'-','color',0.6*[1 1 1]);
        xlabel('t [ d ]');
        ylabel('h_{el} [ m ]');
        %legend('[0 0]','[-4 0]','[-8 0]','avg')
        legend('lake','29','28','27','26','39','38','36','35');

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
print(fig2,str,'-depsc2')

