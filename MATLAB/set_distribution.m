function [ alpha v, tau ] = set_distribution( alpha,v, tau, crack_property )

  % Returns the combined crack distribution and stiff pore distribution 
  
  % alpha             : alpha input is only aspect ratios of the stiff pores
  % v                 : v input is only concentraions of the stiff pores
  % tau               : tau input is only relaxation time constant of the
  %                   : stiff pores
  % crack_property    : Is the properties of the cracks in the system, 1x4 vector.
  % crack_property(1) : Total crack density
  % crack_property(2) : Mean aspect ratio of the crack disttribution 
  % crack_property(3) : delta parameter to use in the betadistr function
  % crack_property(4) : u parameter which sets the end member of the crack
  %                     distribution

  % 09.03.2012
  % Remy Agersborg 
    
  %WARNING: betadistr function is very unstable. It uses a beta propability
  %density function which is very dependent on input parameters
  [alpha_c v_c ]= betadistr(crack_property(3), crack_property(2), crack_property(4), crack_property(1));
  %remove all zero elements in alpha and v
  ind = find(~v_c);
  alpha_c(ind) = [];
  v_c(ind)     = [];
   
  %Add the crack distribution to the rest of the pore types.
  tau_c = ones(size(alpha_c))* tau(1);             % tau: relaxation time 
  alpha = cat(2,alpha_c, alpha);
  v = cat(2,v_c , v);
  tau = cat(2,tau_c,tau);
    
    