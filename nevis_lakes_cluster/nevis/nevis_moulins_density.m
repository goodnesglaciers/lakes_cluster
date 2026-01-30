function [ni_m,sum_m] = nevis_moulins_density(x_m,y_m,gg,oo,aa,ps)
% find indices and locations of moulins nearest coordinates [x_m y_m]
%   also calculates matrix sum_m [n_m-by-nIJ] to describe Voronoi cell for each moulin [each row
%   corresponds to a moulin, having entries 1 for nodes in its cell, 0.5 for nodes on boundary, and 0
%   for nodes outside] [eg area_m = sum_m*(Dx.*Dy) gives the area of each moulin's cell]
% inputs  
%   x_m,y_m   coordinates of moulins  
%   gg     	  grid structre containing node coordinates nx and ny
%   oo        [optional] options
% outputs 
%   n_m       number of moulins [ie. length(x_m)]
%   x_m,y_m   [n_m-by-1] coordinates of moulin grid points
%   ni_m      [n_m-by-1] indices of moulin nodes
%   sum_m     [n_m-by-nIJ] matrix to sum quantities defined on nodes over each moulin's grid cell
%
% 21 August 2014 : taken from hydro_locate_m2
% 25 August 2016 : add in moulin density option LAS

if nargin<4, oo = struct; end
if ~isfield(oo,'random_moulins'), oo.random_moulins = 0; end % pick random moulin coordinates [value should be desired number of moulins]
if ~isfield(oo,'keep_all_moulins'), oo.keep_all_moulins = 1; end % keep all moulins even if on same node
if ~isfield(oo,'move_moulins'), oo.move_moulins = 0; end % move moulins that are outside domain to nearest position inside domain
if ~isfield(oo,'prescribed_catchments'), oo.prescribed_catchments = 0; end % use prescribed shapefile S to give catchments of moulins
if ~isfield(oo,'density_moulins'), oo.density_moulins = 0; end % set moulin number and location from given density

%% pick random moulin coordinates
if oo.random_moulins
    xl = min(min(gg.nx)); xr = max(max(gg.nx)); 
    yb = min(min(gg.ny)); yt = max(max(gg.ny)); 
    n_m = oo.random_moulins;          %number of moulins
    x_m = xl + (xr-xl)*rand(n_m,1);   %x locations
    y_m = yb + (yt-yb)*rand(n_m,1);   %y locations
%     % or
%     x_m = xl + (xr-xl)*rand(n_m,1).^(1/2);   %x locations
%     y_m = yb + (yt-yb)*rand(n_m,1);   %y locations
%     % or 
%     x_m = xl + (xr-xl)*rand(n_m,1).^(1/2);   %x locations
%     f = @(x) 0.0001*x+x.^5;         %x = 0:0.001:1; plot(x,0.5+0.5*f(x-0.5)/f(0.5));
%     y_m = yb + (yt-yb)*( 0.5+0.5*f(rand(n_m,1)-0.5)/f(0.5) );   %dimensionless y locations
end

%% pick random moulin coordinates based on density equation LAS 25 August 2016
if oo.density_moulins
    xl = min(min(gg.nx)); xr = max(max(gg.nx)); 
    yb = min(min(gg.ny)); yt = max(max(gg.ny)); 
    sb = min(aa.s);  st = max(aa.s);           % min and max surface elvation
    
    % Density linear function
    %d_m = oo.density_moulins_slope;            % density slope with elevation (moulins/m2/m elevation)
    %d_b = oo.density_moulins_intercept;        % density at lowest elevation (moulins/m2)
    %p_m = d_b + d_m.*(linspace(sb*ps.z, st*ps.z, 11));% linear decrease in moulin density with elevation (10 bands)
    
    % L. Andrews Thesis Fig. 2.6 Swiss Camp region moulin densities data 2011
    elev_band = (0:100:1600)';       % Andrews Thesis data through 1100 m a.s.l. 
    p_m = [0, 0, 0.19, 0.22, 0.54, 0.48, 0.35, 0.21, 0.19, 0.16, 0.10, 0.07, 0.04, 0.01, 0.005, 0.001]'./(1e6); % moulin density in moulin/m2
    
    x_m = []; y_m = [];
    for i=1:1:(length(elev_band)-1) 
        %elev_band = linspace(sb*ps.z, st*ps.z, 11);         % band widths in S (m)
        [row, col] = find((aa.s.*ps.z) <= elev_band(i+1) & aa.s.*ps.z > elev_band(i)); % index values in band
        area_band = (length(row))*(gg.Dx(1)*ps.x)*(gg.Dy(1)*ps.x);         % find area of band (m2)
        n_m = round(area_band*p_m(i));                       % number of moulins in elevation band
        if n_m > 0
            index_m = round(length(row).*rand(n_m,1)); index_m(index_m==0) = 1;           
            position = row(index_m); % moulins index in aa.s
            x = reshape(gg.nx,gg.nIJ,1); y = reshape(gg.ny,gg.nIJ,1);   % change to x and y column 
            x_m{i} = x(position); y_m{i} = y(position);                 % store
        else 
        end      
    end
    x_m = cat(1,x_m{:}); y_m = cat(1,y_m{:}); % position of moulins across domain
end    
    
%% find indices close to moulins
nx = gg.nx; ny = gg.ny;
if oo.move_moulins, nx = nx(gg.ns); ny = ny(gg.ns); end % restrict to active nodes
n_m = length(x_m);  
ni_m = zeros(n_m,1);
for i_m = 1:n_m
%     ni_m(i_m) = find(nx>=x_m(i_m) & ny>=y_m(i_m),1,'first');
    [~,tmp] = min((reshape(nx,[],1)-x_m(i_m)).^2+(reshape(ny,[],1)-y_m(i_m)).^2); 
    if oo.move_moulins, ni_m(i_m) = gg.ns(tmp); else ni_m(i_m) = tmp; end
end
if ~oo.keep_all_moulins, [ni_m,~,tmp] = unique(ni_m); end; %ni_m = ni_m(tmp); end % [ 26 August 2014: this doesn't look to be actually getting rid of repeated moulins, line below would be better? ]
% if ~oo.keep_all_moulins, ni_m = unique(ni_m,'stable'); end % [ 26 August 2014: changed from above line, option 'stable' doesn't work on older versions of matlab ]
n_m = length(ni_m);
x_m = gg.nx(ni_m); y_m = gg.ny(ni_m);

%% set up catchments for each moulin
sum_m = zeros(length(ni_m),gg.nIJ);
x = reshape(gg.nx,gg.nIJ,1); y = reshape(gg.ny,gg.nIJ,1);
%coordinates of corners that are far away [to prevent any voronoi cells that matter being unbounded]
big = 0.5*(max(x)-min(x) + max(y)-min(y));
big = 0.2*(max(x)-min(x) + max(y)-min(y)); % change 7 Nov
x_far = [max(x)+big; max(x)+big; min(x)-big; min(x)-big];
y_far = [max(y)+big; min(y)-big; max(y)+big; min(y)-big];
%calculate voronoi cells [ last four correspond to far away corners and are irrelevant ]
[v,c] = voronoin([x_m y_m;x_far y_far]);
%find nodes within each moulin cell
for i = 1:length(ni_m),
    [in] = inpolygon(x,y,v(c{i},1),v(c{i},2));
    sum_m(i,:) = in';
end

%% prescribed catchments
if oo.prescribed_catchments
sum_m = zeros(length(ni_m),gg.nIJ);
for i=1:length(ni_m),
    S = oo.catchments;
    xx = S(i).x(1:end-1);
    yy = S(i).y(1:end-1);
    [in] = inpolygon(x,y,xx(~isnan(xx)),yy(~isnan(yy)));
    sum_m(i,:) = in';
end

%split contibution of nodes that are on the boundary of cells
temp = sum(sum_m,1); temp(temp==0) = inf;
sum_m = sum_m*sparse(1:length(temp),1:length(temp),temp.^(-1));

%% plot voronoi cells [uncomment]
% figure;
%     voronoi([x_m; x_far],[y_m; y_far]);
%     axis([min(x) max(x) min(y) max(y)]);
end
