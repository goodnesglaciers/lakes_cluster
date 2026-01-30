function [tt,vv,info] = nevis_timesteps(t_span,vv,aa,pp,gg,oo,pd,ps)
% timestep nevis variables over time t_span
% Inputs:
%   t_span  start and end times and times for saving output
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   tt      struct containing time series of average quantities
%   vv      struct containing final solution variables
%   info    information about last computation [optional]
%
% IJH 14 August 2014 : largely taken from hydro_timesteps_diag
%       1 May 2015 : add possibility for time averaging between timesteps
%       26 October 2015 : added calculation of elastic sheet volume to tt
%       17 Novmber 2015 : changed update boundaries to allow for aa.phi_b
%       to be prescribed externally
% LAS 19 JULY 2016: edited to run on external machines with nohup lines 41 and 147
% LAS 2 DEC 2019: Added sheet and channel flow to be saved at points of interest

if ~isfield(pp,'phi_s'), pp.phi_s = -inf; end

% TIMESTEPPING OPTIONS
if ~isfield(oo,'include_ice'), oo.include_ice = 0; end;         % update ice velocity
if ~isfield(oo,'change_timestep'), oo.change_timestep = 1; end; % change timestep based on previous iterations/convergence
if ~isfield(oo,'dt_max'), oo.dt_max = 1e0; end                  % maximum allowable timestep
if ~isfield(oo,'dt_min'), oo.dt_min = 1e-6; end                 % mimimum allowable timestep
if ~isfield(oo,'dt_factor'), oo.dt_factor = 2; end              % factor to change timestep by
if ~isfield(oo,'small_iter'), oo.small_iter = 2; end            % increase timestep for this number of iterations or less
if ~isfield(oo,'large_iter'), oo.large_iter = 20; end           % decrease timestep for this number of iterations or more 
if ~isfield(oo,'verb'), oo.verb = 0; end                        % verbose screen output
if ~isfield(oo,'adjust_boundaries'), oo.adjust_boundaries = 1; end  % adjust boundary nodes during timestepping
if ~isfield(oo,'dti_boundaries'), oo.dti_boundaries = 4; end    % frequency of timesteps at which to adjust boundaries
% INPUT OPTIONS
if ~isfield(oo,'input_function'), oo.input_function = 0; end    % use pp.input function(t) for moulin inputs
if ~isfield(oo,'runoff_function'), oo.runoff_function = 0; end  % use pp.runoff function(t) for runoff
if ~isfield(oo,'distributed_input'), oo.distributed_input = 0; end  % use distributed input for runoff even in presence of moulins
% SAVING OPTIONS
if ~isfield(oo,'save_timesteps'), oo.save_timesteps = 0; end    % save solution at timesteps in t_span to directory
%if oo.save_timesteps && ~exist([oo.root,oo.fn],'dir'), mkdir([oo.root,oo.fn]); end % make directory if it doesn't already exist
if oo.save_timesteps && ~exist(oo.fn,'dir'), mkdir(oo.fn); end % make directory if it doesn't already exist (nohup syntax edit)
if ~isfield(oo,'save_pts_all'), oo.save_pts_all = 0; end        % save output at particular nodes pts_ni
if ~isfield(oo,'pts_ni'), oo.pts_ni = []; end                   % nodes at which to save output
if ~isfield(oo,'save_phi_av'), oo.save_phi_av = 0; end          % save average pressure between timesteps


% UPDATE BOUNDARIES
if oo.adjust_boundaries && isfield(vv,'nbdy'),         
    gg = nevis_label(gg,vv.nbdy,oo);               % redefine boundary labels
    if ~isfield(aa,'phi_b'), aa.phi_b = max(aa.phi_0,aa.phi_a); end % assign prescribed pressures if phi_b not in aa
    aa.phi = aa.phi_b(gg.nbdy);                    % boundary conditions
end

% TIMESTEPPING
t_start = t_span(1);
t_stop = t_span(end);
ti_save = 1;
t_save = t_span(ti_save); 
ti = 1;
ti_boundaries = 1;
dt1 = oo.dt;    % suggested timestep
t = t_start;
vv.dt = dt1;
vv.t = t;
tt = struct;
if oo.save_phi_av
phi_av = vv.phi; % initial average pressure
phi_min = vv.phi;
phi_max = vv.phi;
end
[vv2,~,~,~,~,~,~,~,~] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % initial expanded variables vv2
while t<t_stop+oo.dt_min

%     if oo.include_ice
%         % [ this option not operational yet ]
%         N = max(0,phi_0-phi);  %[ ignore negative effective pressures ]
%         
%         [ubar,vbar,eta] = ice_force_balance_new2(b,H,A,N,ubar,vbar,eta,tauxx_bdy,tauyy_bdy,tauxy_bdy,pp,gg,oo,opts); 
%         [ub,vb,us,vs,qx,qy,taubx,tauby,tauxxbar,tauyybar,tauxybar,tauxz,tauyz,ui,vi,tau] = ice_velocities_new2(b,H,A,N,ubar,vbar,eta,tauxx_bdy,tauyy_bdy,tauxy_bdy,pp,gg,oo,opts);
%         us(gg.bdy.elbdy,:) = us(gg.bdy.elbdy+1,:); uaa.phi = max(aa.phi_a(gg.nbdy),pp.phi_s);      %i(gg.bdy.elbdy,:) = ui(gg.bdy.elbdy+1,:);  % [ xslope may not be well defined on elbdy : extrapolate values from next points ]
%        %[~,~,wf,~] = ice_velocities_all(b,H,ub,vb,ui,vi,pp,gg,oo,opts); ws = wf(:,1); % vertical surface velocity
% 

%        if oo.ice_update,
%            % update sliding speed [ assumes pd.u was taken as ps.u so that non-dimensional Ub = 1 ]
%            Ub = ( (oo.nmeanx(:,gg.ein)*ub(gg.ein)).^2 + (oo.nmeany(:,gg.fin)*vb(gg.fin)).^2 ).^(1/2);
%        end
%        if oo.melt_update
%            % update melt rate 
%            m = ((pd.G+ ps.tau0*ps.u*( oo.nmeanx(:,gg.ein)*(taubx(gg.ein).*ub(gg.ein))+oo.nmeany(:,gg.fin)*(tauby(gg.fin).*vb(gg.fin)) ))/pd.rho_w/pd.L)/ps.m; 
%        end
%     end
%     end

    %% update input
    aa = nevis_inputs(t,aa,pp,gg,oo);

    %% time series of average quantities
    tt(ti).ti = ti;
    tt(ti).t = t;
    tt(ti).m = sum(aa.m(gg.ns).*gg.Dx(gg.ns).*gg.Dy(gg.ns));  % basal melt, scaled with ps.m*ps.x^2
    tt(ti).E = sum(aa.E(gg.ns).*gg.Dx(gg.ns).*gg.Dy(gg.ns));  % supraglacial input, scaled with ps.m*ps.x^2
    tt(ti).Q_in = vv2.Q_in;                                % inflow, scaled with ps.Q
    tt(ti).Q_out = vv2.Q_out;                              % outflow, scaled with ps.Q
    tt(ti).Q_outQ = vv2.Q_outQ;                            % channel outflow, scaled with ps.Q
    tt(ti).Q_outq = vv2.Q_outq;                            % sheet outflow, scaled with ps.Q
    tt(ti).phi = mean(vv.phi(gg.ns));                      % mean potential, scaled with ps.phi
    tt(ti).N = mean(aa.phi_0(gg.ns)-vv.phi(gg.ns));        % mean effective pressure, scaled with ps.phi
    tt(ti).hs = sum(vv.hs(gg.ns).*gg.Dx(gg.ns).*gg.Dy(gg.ns)); % total cavity sheet volume, scaled with ps.h*ps.x^2
    tt(ti).S = sum(vv.Sx(gg.ein).*(gg.emean(gg.ein,:)*gg.Dx))+sum(vv.Sy(gg.fin).*(gg.fmean(gg.fin,:)*gg.Dy))+sum(vv.Ss(gg.cin).*(gg.cmean(gg.cin,:)*gg.Ds))+sum(vv.Sr(gg.cin).*(gg.cmean(gg.cin,:)*gg.Dr)); % total channel volume, scaled with ps.S*ps.x
    tt(ti).he = sum(vv2.he(gg.ns).*gg.Dx(gg.ns).*gg.Dy(gg.ns)); % total elastic sheet volume, scaled with ps.h*ps.x^2
    if oo.save_pts_all
         tt(ti).pts_phi = vv.phi(oo.pts_ni);
         tt(ti).pts_N = aa.phi_0(oo.pts_ni)-vv.phi(oo.pts_ni);
         tt(ti).pts_hs = vv.hs(oo.pts_ni);
         tt(ti).pts_he = vv2.he(oo.pts_ni);
         tt(ti).pts_Sx = vv.Sx(oo.pts_ni);
         tt(ti).pts_Sy = vv.Sy(oo.pts_ni);
         tt(ti).pts_Ss = vv.Ss(oo.pts_ni);
         tt(ti).pts_Sr = vv.Sr(oo.pts_ni);
%          tt(ti).pts_qQ = vv.qQ(oo.pts_ni);    % sheet discharge, scaled with ps.Q % LAS 2 DEC 2019
%          tt(ti).pts_Q = vv.Q(oo.pts_ni);    % channel discharge, scaled with ps.Q % LAS 2 DEC 2019
         if oo.include_ice
         tt(ti).pts_us = vv.us(oo.pts_ni);
         end
    end
    
    %% saving    
    if t >= t_save,
        if oo.save_phi_av,
            vv.phi_av = phi_av;
            vv.phi_min = phi_min;
            vv.phi_max = phi_max;
            phi_av = vv.phi; % reset average pressure
            phi_min = vv.phi;
            phi_max = vv.phi;
        end
        vv.ti = ti;
        vv.ti_save = ti_save;
        fn = [oo.root,oo.fn];
        disp(['nevis_timesteps: Saving ',fn,' ...']);
        save(fn,'vv','tt','-append');
        disp('nevis_timesteps: Done');
        % save full solution
        if oo.save_timesteps
            %fnt = [oo.root,oo.fn,'/',int2four(ti_save)]; % old line
	    fnt = [oo.fn,'/',int2four(ti_save)]; % LAS nohup edit	
            disp(['nevis_timesteps: Saving ',fnt,' ...']);
            save(fnt,'vv','tt','pp','pd','ps'); 
            disp('nevis_timesteps: Done');
        end
        % [ add other saving options here ]
        if t>t_stop-oo.dt_min, break; end;
        ti_save = ti_save+1;
        t_save = t_span(ti_save);  
    elseif t > t_stop-oo.dt_min, break; 
    end;

    %% timestep hydrology
    tic;
    accept = 0;
    decreased = 0;
    increased = 0;
    while ~accept
        dt = dt1; %suggested timestep
        if t + dt > t_save-oo.dt_min, dt = t_save-t; end % adjust timestep if this would take too close to t_save
        disp(['nevis_timesteps: dt = ',num2str(dt),' ...']);
        [vv1,vv2,info] = nevis_timestep(dt,vv,aa,pp,gg,oo);

        %% check success and adjust size of timestep [ taken from hydro_timestep_diag ]
        if ~info.failed, 
            accept = 1; 
            t = t + dt;
            vv = vv1;
            vv.dt = dt;
            vv.t = t;
            comp_time = toc;
            disp(['nevis_timesteps: Done [ ',num2str(comp_time),' s, ',num2str(info.iter_new-1),' iterations ]']);
        end
        if oo.change_timestep && info.iter_new-1 <= oo.small_iter && ~info.failed && dt1 < oo.dt_max && ~decreased
            dt1 = min(oo.dt_factor*dt1,oo.dt_max); increased = 1;
            if oo.verb, disp(['nevis_timesteps: Increase suggested timestep to ',num2str(dt1)]); end
            continue;
        elseif oo.change_timestep && ( info.iter_new-1 >= oo.large_iter || info.failed ) && dt1 > oo.dt_min && ~increased
            dt1 = max(dt1/oo.dt_factor,oo.dt_min); decreased = 1;
            if oo.verb, disp(['nevis_timesteps: Decrease suggested timestep to ',num2str(dt1)]); end
            continue;
        end
        if ~accept
            comp_time = toc;
            disp(['nevis_timesteps: Failed [ ',num2str(comp_time),' s, ',num2str(info.iter_new-1),' iterations, dt = ',num2str(dt),' ]']);
            return;
        end

    end
    
    %% time averaging of pressure
    if oo.save_phi_av,
        phi_av = ( phi_av*(t-t_span(ti_save-1)-dt)+vv.phi*dt )/(t-t_span(ti_save-1));
        phi_min = min(phi_min,vv.phi);
        phi_max = max(phi_max,vv.phi);
    end
    
    %% check and adjust boundary nodes
    %oo.verb=1;
    if oo.adjust_boundaries && ti>=ti_boundaries
    ti_boundaries = ti+oo.dti_boundaries;
    ni1 = gg.nbdy(vv2.R_bdy<0); % Dirichlet nodes with inflow
    if ~isfield(gg,'n1m'), gg.n1m = gg.n1; end % boundary nodes adjacent to margin
    ni2 = gg.n1m(vv.phi(gg.n1m)-aa.phi_b(gg.n1m)>pp.p_a_reg); % boundary nodes with too high pressure (LAS 18 nov 2015 phi_0)
    if ~isempty(ni1) || ~isempty(ni2)
        if ~isempty(ni1), 
            if oo.verb, disp('nevis_timesteps: Removing Dirichlet indices ...');
            disp(ni1); end
        end
        if ~isempty(ni2),
            if oo.verb, disp('nevis_timesteps: Adding Dirichlet indices ...');
            disp(ni2); end
        end
        vv.nbdy = union(setdiff(gg.nbdy,ni1),ni2);
        gg = nevis_label(gg,vv.nbdy,oo);                        % redefine boundary labels
        % aa.phi = max(aa.phi_a(gg.nbdy),pp.phi_s);             % boundary conditions
        aa.phi = aa.phi_b(gg.nbdy);                             % boundary conditions for phi_b
        % nevis_plot_grid(gg,gg.nbdy);  % for maintenance
    end
    end 
    
    disp(['nevis_timesteps: t = ',num2str(t),' [ / ',num2str(t_stop),' ]']);
    ti = ti + 1; 
    
end
end
