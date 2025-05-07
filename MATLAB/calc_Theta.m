function [Theta ] = calc_Theta( v, omega, gamma,tau, Kd_uuvv, Dr,k, kappa, kappa_f)

% Returns the Theta tensor (6x6 matrix) for more explanation see e.g.
% Agersborg et al. 2009 or "The effects of drained and undrained loading in
% visco-elsatic waves in rock-like composites" M. Jakobsen and T.A. % Johansen. 
%(2005). Int. J. Solids and Structures (42). p. 1597-1611

% v       : concentration of all the empty cavities (1x(numbers of empty cavities)
%           vector) 
% omega   : frequency (2*pi*f)
% gamma   : gamma factor of all the inclusions (1x(numbers of empty cavities)
%           vector) 
% tau     : Relaxation time constant (1x(numbers of empty cavities)
%           vector)
% Kd_uuvv : Tensor sum of the dry K tensor for all the cavities (1x(numbers of empty cavities)
%           vector) 
% Dr      : Permeability/viscosity
% k       : Wave number vector
% kappa   : Bulk modulus of the host material.
% kappa_f : Bulk modulus of the fluid.


% 09.03.2012
% Remy Agersborg 

sigma_a = sum( (v./(1+i*omega.*gamma.*tau)));
sigma_b = sum( (v./(1+i*omega.*gamma.*tau)).*Kd_uuvv );

Theta = kappa_f/((1-kappa_f/kappa)*sigma_a + kappa_f*sigma_b - (i*k*k/omega)*Dr*kappa_f);
