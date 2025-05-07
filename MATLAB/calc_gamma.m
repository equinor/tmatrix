function [ gamma ] = calc_gamma(kappa_f,kappa, Kd_uuvv )
% Returns the gamma factor for inclusion.

% kappa_f : Bulk modulus of the fluid.
% kappa   : Bulk modulus of the host material.
% Kd_uuvv : Sum of the dry K-tensor for inclusion under consideration.

% Equations used can be found in:
% Agersborg (2007), phd thesis:
% https://bora.uib.no/handle/1956/2422


% 09.03.2012
% Remy Agersborg 

gamma = 1 - kappa_f/kappa + kappa_f.*Kd_uuvv;
end


