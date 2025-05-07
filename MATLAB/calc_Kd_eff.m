function [ Kd_eff_connected, Kd_eff_isolated ] = calc_Kd_eff(C0,kappa_f, alpha, v, Gd, ctrl )
% Returns the effective dry K-tensor (6x6x(numbers of inclusions) matrix.
% If there is no connected or no isolated pores, the function returns a NaN for
% the case which is not considered. E.g. if only isolated pores, the Kd_eff_connected = NaN. 

% C0      : Stiffness tensor of the host material (6x6 matrix).
% kappa_f : bulk modulus of the fluid.
% alpha   : aspect ratio  (structur: alpha.isolated and alpha.connected).
% v       : Concentration (structur: v.isolated and v.connected).
% Gd      : correlation function (6x6 matrix).
% ctrl    : Control parameter (vector).
% ctrl(2) = 0 :only isolated pores.
% ctrl(2) = 1 :both isolated and connected pores.
% ctrl(2) = 2 :only connected pores.

%Note: When isolated pores, the pores are considered as filled when
%calculating the dry effective K-tensor. 

% Equations used can be found in:
% Agersborg (2007), phd thesis:
% https://bora.uib.no/handle/1956/2422

% 09.03.2012
% Remy Agersborg 

C1_isolated = zeros(6,6);
C1_connected = zeros(6,6);
if ctrl(2) ~= 2
    C1_isolated  = calc_isolated_part(C0,kappa_f,alpha.isolated, v.isolated, ctrl);
end
if ctrl(2) ~= 0
    C1_connected = calc_isolated_part(C0,0,alpha.connected, v.connected, ctrl);
end

C1dry = C1_isolated + C1_connected;

C2dry = C1_connected*Gd*C1_connected + C1_connected*Gd*C1_isolated + C1_isolated*Gd*C1_connected + C1_isolated*Gd*C1_isolated;
 
I4 = eye(6);
C_eff_dry = C0 + C1dry/(I4 + (C1dry)\C2dry);

temp = I4/((I4 + (C1dry)\C2dry) )/(C_eff_dry);

%if only connected or mixed connected and isolated
if ctrl(2) ~= 0
    Kd_eff_connected = zeros(6,6,length(alpha.connected));
    for j = 1: length(alpha.connected)
        G = Gtensor(C0, alpha.connected(j));
        Kd_eff_connected(:,:,j) = ((I4 +G*C0))\temp ;
    end
end
if ctrl(2) ~= 2
    Kd_eff_isolated = zeros(6,6,length(alpha.isolated));
    for j = 1: length(alpha.isolated)
        G = Gtensor(C0, alpha.isolated(j));
        Kd_eff_isolated(:,:,j) = ((I4 +G*C0))\temp;
    end
end   

if ctrl(2) == 2
    Kd_eff_isolated = NaN;
end
if ctrl(2) == 0
    Kd_eff_connected = NaN;
end
    

end
