function N = nevis_uplift_pressure(h,pp)
% finds effective pressure N corresponding to elastic uplift h, using
%   h = c_e_reg2*(-N+N_reg2/2) for N<0
%   h = c_e_reg2*N_reg2/2*(1-N/N_reg2)_+.^2 for N>0
% those points with h = 0 have N>N_reg2 - they are assigned NaN
% Inputs: 
%    h elastic sheet uplift
%    pp parameter struct
% Ouputs: 
%    N effective pressure corresponding to h

% if non-dimensional
c_e_reg2 = pp.c31;  % scaled with ps.phi/ps.h
N_reg2 = pp.c32; % sclaed with ps.phi

% % if dimensional
% c_e_reg2 = pd.c_e_reg2;
% N_reg2 = pd.N_reg2;


N = N_reg2.*max(0,1-h/(c_e_reg2*N_reg2/2)).^(1/2); % effective pressure in regularization range 0<N<N_reg2
ii = (h>c_e_reg2*N_reg2/2); % overburden nodes
N(ii) = N_reg2/2-h(ii)/c_e_reg2; 
N(h<=0) = NaN; % unsconstrained pressure N>N_reg2

end

