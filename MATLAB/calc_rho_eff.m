function [rho_eff] = calc_rho_eff(rho_f, rho_s, phi)
% Returns the effective density 
% 
% rho_f : density of the fluid
% rho_s : density of the solid
% phi   : porosity

% 09.03.2012
% Remy Agersborg 

rho_eff = phi*rho_f + (1-phi)*rho_s;

end