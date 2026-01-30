function [x_tran,y_tran,s_tran,var_tran] = nevis_var_transect_xe(x,y,N,var,gg)
% [x_tran,y_tran,phi_tran] = nevis_phi_transect(x1,y1,x2,y2,N,phi,gg)
%
% Interpolate pressure along a straight line segment between coordinates x
% and y, with N points
% Inputs:
%   x,y           coordinates of end points of transect (2-by-1 vectors, eg
%   obtained from nevis_get_transect)
%   N             number of points to include on transect
%   var           variable of interest
% Outputs:
%   x_tran,y_tran coordinates of transect
%   s_tran        arc length along transect
%   var_tran      pressure on transect (interpolated)
%   
% 02 Nov 2015 LAS

% reshape pressure onto a matrix
var = reshape(var,gg.eI,gg.eJ);  
xx = reshape(gg.ex,gg.eI,gg.eJ);
yy = reshape(gg.ey,gg.eI,gg.eJ);

% intepolate
x_tran = linspace(x(1),x(2),N);
y_tran = linspace(y(1),y(2),N);
dx = x_tran(2:end)-x_tran(1:end-1);
dy = y_tran(2:end)-y_tran(1:end-1);
s_tran = [0 cumsum((dx.^2+dy.^2).^(1/2))];
var_tran = interp2(xx',yy',var',x_tran,y_tran);

end