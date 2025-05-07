function [ alpha_m, v_m, tau_m, gamma_m ] = calc_pressure( alpha, v, C0, Gd, dP, tau, gamma, kappa_f, ctrl)
% Returns the altered aspect ratios and corresponding concentrations due to change of effective pressure dP . If any
% inclusions are closed due to the pressure change, the length of tau and gamma is also
% changed to correspond to the new length of alpha and v. 

% NOTE: THIS IS ONLY VALID FOR DRAINED CONDITION 
 
% alpha   : aspect ratios of the inclusions (structure: alpha.isolated and
%           alpha.connected)
% v       : concentration of the inclusions (structure: v.isolated and
%           v.connected)
% C0      : stiffness tensor of the host material
% Gd      : the correlation function (green's tensor 6x6 matrix) .
% dP      : change in effective pressure.
% tau     : relaxation time constant ((1x(numbers of connected pores)
%           vector)
% gamma   : gamma factor ((1x(numbers of connected pores) vector)
% kappa_f : Bulk modulus of the fluid
% ctrl    : Control parameter 
% ctrl(2) = 0 :only isolated pores
% ctrl(2) = 1 :both isolated and connected pores
% ctrl(2) = 2 :only connected pores

% Equations used can be found in:
% Agersborg (2007), phd thesis:
% https://bora.uib.no/handle/1956/2422

% 09.03.2012
% Remy Agersborg 

alpha_threshold = 0;
v_threshold     = 0;
% find the effective dry K-tensors for the system
[ Kd_eff_connected, Kd_eff_isolated ] = calc_Kd_eff(C0, kappa_f, alpha, v, Gd, ctrl );
if ctrl(2) ~= 2
    count_isolated   = size(Kd_eff_isolated);
    v_n_isolated     = zeros(1, count_isolated(3));
    alpha_n_isolated = zeros(1, count_isolated(3));
end
if ctrl(2) ~= 0
    count_connected   = size(Kd_eff_connected);
    v_n_connected     = zeros(1, count_connected(3));
    alpha_n_connected = zeros(1, count_connected(3));
end
sum = 0.d0;

%find the sum in the eq. 21 Jakobsen and Johansen 2005.
if ctrl(2) ~= 2
    for j=1: count_isolated(3)
        sum = sum + v.isolated(j)*(Kd_eff_isolated(1,1,j)+Kd_eff_isolated(1,2,j)+Kd_eff_isolated(1,3,j)  ...
                                  +Kd_eff_isolated(2,1,j)+Kd_eff_isolated(2,2,j)+Kd_eff_isolated(2,3,j)  ...
                                  +Kd_eff_isolated(3,1,j)+Kd_eff_isolated(3,2,j)+Kd_eff_isolated(3,3,j));
    end;
end
if ctrl(2) ~= 0
    for j=1: count_connected(3)
        sum = sum + v.connected(j)*(Kd_eff_connected(1,1,j)+Kd_eff_connected(1,2,j)+Kd_eff_connected(1,3,j)  ...
                                  +Kd_eff_connected(2,1,j)+Kd_eff_connected(2,2,j)+Kd_eff_connected(2,3,j)  ...
                                  +Kd_eff_connected(3,1,j)+Kd_eff_connected(3,2,j)+Kd_eff_connected(3,3,j));
    end;
end
%Find the new concentration of inclusion
if ctrl(2) ~= 2
    for j = 1:count_isolated(3)
        v_new = v.isolated(j)*(1 -((Kd_eff_isolated(1,1,j)+Kd_eff_isolated(1,2,j)+Kd_eff_isolated(1,3,j) ...
                    +Kd_eff_isolated(2,1,j)+Kd_eff_isolated(2,2,j)+Kd_eff_isolated(2,3,j) ...
                    +Kd_eff_isolated(3,1,j)+Kd_eff_isolated(3,2,j)+Kd_eff_isolated(3,3,j)) - sum)*dP) ;
        if v_new <= v_threshold 
            v_n_isolated(j) = nan;
        else
            v_n_isolated(j) = v_new;
        end
    end
end
if ctrl(2) ~= 0
    for j = 1:count_connected(3)
        v_new = v.connected(j)*(1 -((Kd_eff_connected(1,1,j)+Kd_eff_connected(1,2,j)+Kd_eff_connected(1,3,j) ...
                    +Kd_eff_connected(2,1,j)+Kd_eff_connected(2,2,j)+Kd_eff_connected(2,3,j) ...
                    +Kd_eff_connected(3,1,j)+Kd_eff_connected(3,2,j)+Kd_eff_connected(3,3,j)) - sum)*dP) ;
        if v_new <= v_threshold 
            v_n_connected(j) = nan;
        else
            v_n_connected(j) = v_new;
        end
    end
end

%Find the new alpha of inclusion
if ctrl(2) ~= 2
    for j = 1:count_isolated(3)
        alpha_new = alpha.isolated(j)*(1 -(Kd_eff_isolated(3,1,j)+Kd_eff_isolated(3,2,j)+Kd_eff_isolated(3,3,j) ...
                                        - Kd_eff_isolated(1,1,j)-Kd_eff_isolated(1,2,j)-Kd_eff_isolated(1,3,j))*dP) ;
        if alpha_new <= alpha_threshold  
            alpha_n_isolated(j) = nan;
        else
            alpha_n_isolated(j) = alpha_new;
        end
    end
end
if ctrl(2) ~= 0
    for j = 1:count_connected(3)
        alpha_new = alpha.connected(j)*(1 -(Kd_eff_connected(3,1,j)+Kd_eff_connected(3,2,j)+Kd_eff_connected(3,3,j) ...
                                        - Kd_eff_connected(1,1,j)-Kd_eff_connected(1,2,j)-Kd_eff_connected(1,3,j))*dP) ;
        if alpha_new <= alpha_threshold
            alpha_n_connected(j) = nan;
        else
            alpha_n_connected(j) = alpha_new;
        end
    end
end

%Check for nan enteries which will be removed
if ctrl(2) ~= 2
    [alpha_m.isolated, v_m.isolated ] =  check_for_nan_enteries_isolated(alpha_n_isolated ,v_n_isolated);
end
if ctrl(2) ~= 0
   [alpha_m.connected, v_m.connected, tau_m, gamma_m ] =  check_for_nan_enteries_connected(alpha_n_connected ,v_n_connected, tau, gamma);
end



end




    
    




