% voronoi plots for moulin catchments

% location of each moulin
x_m = gg.nx(pp.ni_m); y_m = gg.ny(pp.ni_m);

% catchments for each moulin
sum_m = zeros(length(pp.ni_m),gg.nIJ);
x = reshape(gg.nx,gg.nIJ,1); y = reshape(gg.ny,gg.nIJ,1);

% voronoi cells
big = 0.2*(max(x)-min(x) + max(y)-min(y));
x_far = [max(x)+ big; max(x) + big; min(x)-big; min(x)-big];
y_far = [max(y)+big; min(y)-big; max(y)+big; min(y)-big];
%calculate voronoi cells [ last four correspond to far away corners and are irrelevant ]
[v,c] = voronoin([x_m y_m;x_far y_far]);

%find nodes within each moulin cell
%for i = 1:length(pp.ni_m),
%    [in] = inpolygon(x,y,v(c{i},1),v(c{i},2));
%    sum_m(i,:) = in';
%end
radius=6378137.0; eccen=0.08181919; lat_true=70; lon_posy=-45; % projection parameters
foot = [68.884795, -50.182106; 68.952545, -49.409098;
        %68.393966, -48.867272; 68.282561, -49.831050];
        68.469744, -48.970285; 68.399076, -49.709055];
moulin = [68.723585, -49.536195]; % moulin is (0,0)    
[moulin_x,moulin_y] = polarstereo_fwd(moulin(1),moulin(2),radius,eccen,lat_true,lon_posy);
[foot_x,foot_y] = polarstereo_fwd(foot(:,1),foot(:,2),radius,eccen,lat_true,lon_posy);
footprint(:,1) = foot_x-moulin_x; footprint(:,2) = foot_y - moulin_y;

FLXX = [68.7735, 310.0752-360;
        68.7559, 310.2587-360;
        68.7373, 310.4469-360;
        68.7253, 310.5542-360;
        68.7158, 310.8561-360;
        68.7192, 311.1320-360];
[flxx_x,flxx_y] = polarstereo_fwd(FLXX(:,1),FLXX(:,2),radius,eccen,lat_true,lon_posy);
flxx_plot(:,1) = flxx_x-moulin_x; flxx_plot(:,2) = flxx_y - moulin_y;


figure;
     voronoi((ps.x/10^3).*[x_m; x_far],(ps.x/10^3).*[y_m; y_far],'k');
     hold on
    
    % add moulins
    if ~isfield(pp,'ni_m'), pp.ni_m = []; end
    mscale = 100;  %amount to scale moulin input by to make size of moulin dot
    x = (ps.x/10^3)*gg.nx(pp.ni_m);
    y = (ps.x/10^3)*gg.ny(pp.ni_m);
    hold on;
    for i_m = 1:length(pp.ni_m),
        if aa.E(pp.ni_m(i_m))>0,
            plot(x(i_m),y(i_m),'ko','Markersize',4+aa.E(pp.ni_m(i_m))/mscale,'MarkerFaceColor',1*[1 1 1]); % mark moulins   
        else
            plot(x(i_m),y(i_m),'ko','Markersize',4,'MarkerFaceColor',0.8*[1 1 1]); % mark moulins  
        end
    end
    
    % add outlets
    plot((ps.x/10^3)*gg.nx(gg.nbdy),(ps.x/10^3)*gg.ny(gg.nbdy),'ko','markersize',4,'markerfacecolor','y');
    
    % add boundary
    % plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5); 
    
    % add Joughin 2013 footprint
    patch(footprint(:,1)./1000, footprint(:,2)./1000,'k','FaceColor','none'); 
   axis([-80 60 -50 30]);