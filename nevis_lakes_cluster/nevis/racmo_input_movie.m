% 06 JULY 2016
% RACMO DISTRIBUTED INPUT CHECK VIDEO

format compact;
clear
oo.root = '';           % filename root
oo.fn = mfilename;      % filename
oo.code = './nevis';   % code directory
addpath(oo.code);

%% parameters
% default parameters
[pd,oo] = nevis_defaults([],oo);    
% alter default parmaeters 
pd.c_e_reg2 = 0.01/1e3/9.81;        % elastic sheet thickness [m/Pa]
pd.N_reg2 = 1e4; % 1e3              % regularisation pressure for elastic sheet thickness 
pd.u_b = 100/pd.ty;                 % sliding speed [m/s]
pd.sigma = 1e-3;                    % englacial void fraction
pd.h_r = 0.1; % .1                  % roughness height [m]
pd.l_r = 10; % 10                   % roughness length [m]
pd.k_s = 1e-4;                      % sheet permeability constant
pd.melt = pd.G/pd.rho_w/pd.L;       % geothermal heat derived basal melt [m/s]
% non-dimensionalise
[ps,pp] = nevis_nondimension(pd);   

%% grid and geometry
load('morlighem_for_nevis_200km'); % load bedmap
dd = morlighem_for_nevis_200km;
dd.skip = 6;
gg = nevis_grid(dd.X_km(1:dd.skip:end,1)/ps.x,dd.Y_km(1,1:dd.skip:end)/ps.x,oo); 
b = reshape(dd.B_km(1:dd.skip:end,1:dd.skip:end)/ps.z,gg.nIJ,1);
s = reshape(dd.S_km(1:dd.skip:end,1:dd.skip:end)/ps.z,gg.nIJ,1);
load('runoff_2009_nevis200.mat'); % load data (previously collated)

%% animate 
% Preallocate movie structure.
%mov(1:36) = struct('cdata', [],'colormap', parula);

 for i=1:1:365
    figure(1)
    clf
    contourf((reshape(runoff_2009_nevis200(i,:),gg.nI,gg.nJ))',(0:5:50))
    colorbar
    caxis([0 50])
    title(num2str(i))
    
    drawnow
    F(i) = getframe(gcf);
 end

fig = figure;
movie(fig,F,1);
%movie2avi(mov,'racmo2009_nevis200km.avi')