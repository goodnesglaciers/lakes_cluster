%[Nnev] = prep_nevis_for_coherence(fn, DOI);
%% Pepare nevis output values for mtm analysis
%% Laura A. Stevens 9 Mar 2017
%%  EDIT 14 Mar 2017 -- 1 km grid
close all; 

%% Load nevis and terrasarx files of interest
load nevis_170207b.mat
load zz_nevis.mat
load Nwinter.mat
load('../../Joughin_Cryosphere2013/terravel.mat')

%% nevis reshaping
%get rid of points outside domain
nx(gg.nout) = NaN; ex(gg.eout) = NaN; fx(gg.fout) = NaN; cx(gg.cout) = NaN;
% get x and y fields
xx = (ps.x/10^3)*gg.nx; xxx = reshape(xx, size(xx,1)*size(xx,2), 1);
yy = (ps.x/10^3)*gg.ny; yyy = reshape(yy, size(yy,1)*size(yy,2), 1);

% interpolate nevis to terrasarX grid
terraxx = (terravel.winter.vel.*(terravel.polar.xx_moulin./1000))./(terravel.winter.vel);
terrayy = (terravel.winter.vel.*(terravel.polar.yy_moulin./1000))./(terravel.winter.vel);
F = scatteredInterpolant(xxx,yyy,reshape(Nwinter,size(xx,1)*size(xx,2), 1)); 
N_ter(:,:,1) = F(terraxx, terrayy);
for i=2:9
    F = scatteredInterpolant(xxx,yyy,reshape(zz_nevis(i-1).N,size(xx,1)*size(xx,2), 1)); 
    N_ter(:,:,i) = F(terraxx, terrayy);
end

%% rotate into a size that works for everyone...MUST BE SQUARE, No NaNs
% shift couterclockwise wise in degrees
theta = -10.4*(pi/180); 
R = [cos(theta) -sin(theta); sin(theta) cos(theta)]; % rotation matrix
v = [reshape(terraxx, size(terraxx,1)*size(terraxx,2), 1)';...
    reshape(terrayy, size(terraxx,1)*size(terraxx,2), 1)'];
so = (R*v)'; 
SO = so; SO(any(isnan(SO),2),:) = []; % remove NaNs
SO = sortrows(so,1);  % sort based on x value
% create new terrasarx X and Y meshgrid from the rotated positions at 1 km
% resolution
[txx, tyy] = meshgrid(min(SO(:,1))+1:1:max(SO(:,1)), min(SO(:,2)):1:max(SO(:,2))); 

% square matrix xy
width = (size(txx, 2));
step = 1;
txxb = txx((1:step:width),(1:step:width)); tyyb = tyy((1:step:width),(1:step:width));
txxb = txxb -  min(min(txxb)); % put 0 at 
tyyb = tyyb - min(min(tyyb));

for i=1:9
    all = horzcat(so, reshape(N_ter(:,:,i), size(terraxx,1)*size(terraxx,2), 1)); % rotated
    all(any(isnan(all),2),:) = []; % remove NaNs
    F = scatteredInterpolant(all(:,1),all(:,2),all(:,3));
    
    % rotated N
    N_rot(:,:,i) = F(txx, tyy);
    % rotated N into a square matrix
    N_rot_sqr(:,:,i) = N_rot((1:step:width),(1:step:width),i);
    % as a difference from winter N
    N_rot_sqr_difference(:,:,i) = N_rot_sqr(:,:,i)-N_rot_sqr(:,:,1); % current less winter
    
    % rotated N into a square matrix
    Nup_rot_sqr(:,:,i) = N_rot((1:step:width),(width:step:end),i);
    % as a difference from winter N
    Nup_rot_sqr_difference(:,:,i) = Nup_rot_sqr(:,:,i)-Nup_rot_sqr(:,:,1); % current less winter
    
end

%% save data
Nnev.N_rot = N_rot;
Nnev.N_rot_sqr = N_rot_sqr;
Nnev.N_rot_sqr_difference = N_rot_sqr_difference;
Nnev.Nup_rot_sqr = Nup_rot_sqr;
Nnev.Nup_rot_sqr_difference = Nup_rot_sqr_difference;
Nnev.txx = txx;
Nnev.txxb = txxb;
Nnev.tyy = tyy;
Nnev.tyyb = tyyb;
save Nnev_1km.mat Nnev
