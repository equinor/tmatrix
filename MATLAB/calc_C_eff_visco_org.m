function [C_eff ] = calc_C_eff_visco(Vs, K_r, eta_f, v, gamma, tau, Kd_uuvv, kappa, kappa_f, C0,C1, td_bar, X_bar, Gd, frequency)
%    Returns the effective stiffness tensor C* for a visco-elastic system (6x6xnumber of frequencies).

%    Vs        : The velocity used to calculate the wave number 
%    K_r       : Klinkenberg permability
%    eta_f     : viscosity (P)
%    v         : concentration of the inclusions which are connected with
%                respect to fluid flow
%    gamma     : gamma factor for each inclusion (1x(number of connected
%                inclusions) vector)
%    tau       : relaxation time constant
%    Kd_uuvv   : Kd_uuvv for each connected inclusion (1x(number of connected
%                inclusions) vector
%    kappa     : Bulk modulus of host material
%    kappa_f   : Bulk modulus of the fluid
%    C0        : The stiffness tensor of host material (6x6 matrix)
%    C1        : First order correction matrix(6x6 matrix). If there are
%                isolated inclusions, C1 is sum of concentration and t-matrices of the isolated part of the porosity 
%    td_bar    : t-matrices of the connected inclusions(6x6x(numbers of
%                inclusions) matrix)
%    X_bar     : X-tensor of the connected inclusions (6x6x(numbers of inclusions) matrix)
%    Gd        : correlation function (6x6 matrix).
%    frequency : frequency under consideration (1x(number of frequencies)
%                vector)

% Equations used can be found in:
% Agersborg (2007), phd thesis:
% https://bora.uib.no/handle/1956/2422


% 09.03.2012
% Remy Agersborg 

    Dr = K_r/eta_f;

    C_eff = zeros(6, 6, length(frequency));
    
    for j=1:length(frequency)
        omega = 2*pi*frequency(j);
        k = omega/Vs;
        Theta = calc_Theta(v, omega, gamma,tau, Kd_uuvv, Dr,k, kappa, kappa_f);
        Z_bar = calc_Z(C0,td_bar, omega, gamma, v, tau );
        t_bar = calc_t(td_bar, Theta, X_bar, Z_bar, omega, gamma,tau,kappa_f);
       
        C_eff(:,:,j) = calc_C_eff(C0, C1, Gd, t_bar, v);
    end


end






