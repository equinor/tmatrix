function [ alpha_final, v_final, tau_final ] = sort_pore_types( alpha, v, tau, per_inc_ani )
% If there is both anisotropic and isotropic pore types in the system (0 < per_inc_ani < 1)this
% function divide alpha (aspect ratios) and v (concentration)into a
% isotropic and anisotropic part.


alpha_final = cat(2,alpha,alpha);
v_final     = cat(2,(1-per_inc_ani)*v, per_inc_ani*v);
tau_final   = cat(2,tau, tau);

end