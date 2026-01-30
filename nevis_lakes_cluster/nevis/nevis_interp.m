function out = nevis_interp(in,x_in,y_in,proj_in,x,y,proj,oo)
% interpolate quantity in defined on mesh x_in,y_in in projection proj_in, to mesh x,y in projection
% proj
%
% Example interpolate speed v defined on grid [x1,y1] in proj1 onto grid [x,y] in proj2
%  proj1.type = 'joughin'; proj2.type = 'bamber';
%  v = nevis_interp(v,x_in,y_in,proj1,x,y,proj2);
%
% 23 Sept 2014: modified from Mauro Werder's interp_geo.m [ might want to include option for
%   temporal interpolation here too ]

if nargin<8, oo = struct; end
if ~isfield(oo,'structured_input'), oo.structured_input = 0; end
if size(x_in)~=size(in), [x_in,y_in] = meshgrid(x_in,y_in); end % meshgrid input if x_in and y_in included as vectors?

% transformation to input grid
[x2,y2] = nevis_transform(x,y,proj,proj_in);

if oo.structured_input
% structured input grid
out = interp2(x_in,y_in,in,x2,y2,'linear',NaN);
else
% unstructured input grid
out = griddata(x_in,y_in,in,x2,y2,'linear');
end

end