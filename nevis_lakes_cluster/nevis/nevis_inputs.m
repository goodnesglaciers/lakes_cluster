function aa = nevis_inputs(t,aa,pp,gg,oo)
% aa = nevis_inputs(t,aa,pp,gg,oo)
% update inputs in aa for time t
%
% 21 August 2014
% 11 Dec 2014 : note that if runoff_function includes interp2, it may not work on earlier matlab 
%   versions if interpolating NaN values
% 12 November 2015 : edited to use oo.distributed_input option whether
%    oo.input_function is used or not (rather than only if it's not) by
%    reordering
% 11 May 2016 : edited to not use oo.distributed_input option when
%    oo.runoff_function is used (LAS)

if nargin<5, oo = struct; end
if ~isfield(oo,'input_function'), oo.input_function = 0; end    % use pp.input function(t) for moulin inputs
if ~isfield(oo,'runoff_function'), oo.runoff_function = 0; end  % use pp.runoff function(t) for runoff
if ~isfield(oo,'distributed_input'), oo.distributed_input = 0; end  % use distributed input for runoff even in presence of moulins

    
    if oo.runoff_function
        r = pp.runoff_function(t); r(gg.nout) = 0;
    else
        r = (1-pp.E_amp*cos(2*pi*t/pp.td)).*max(pp.meltE(t)-pp.E_lapse*aa.s,0); r(gg.nout) = 0;
    end
%% old lines    
%     if isfield(pp,'ni_m') && ~isempty(pp.ni_m) && ~oo.distributed_input
%         aa.E = 0*ones(gg.nIJ,1); 
%         aa.E(pp.ni_m) = (pp.sum_m*(r.*gg.Dx.*gg.Dy))./gg.Dx(pp.ni_m)./gg.Dy(pp.ni_m);
%     else
%         aa.E = r; 
%     end

% %% 11 May 2016 edits. Use distributed input only when no runoff input
    if isfield(pp,'ni_m') && ~isempty(pp.ni_m) && oo.distributed_input && ~oo.runoff_function
        aa.E = 0*ones(gg.nIJ,1); 
        aa.E(pp.ni_m) = (pp.sum_m*(r.*gg.Dx.*gg.Dy))./gg.Dx(pp.ni_m)./gg.Dy(pp.ni_m);
    else
        aa.E = r; 
    end
%%    
    if oo.input_function
        aa.E(pp.ni_m) = (pp.input_function(t))./gg.Dx(pp.ni_m)./gg.Dy(pp.ni_m);
    else

    end
    % boundary input
    if isfield(pp,'Q_in')
        aa.Qx = 0*ones(length(gg.ebdy),1);
        aa.Qx(pp.ei_in) = pp.Q_in(t);
    end
end