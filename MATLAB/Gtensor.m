function G = Gtensor(C0, alpha)

%Returns the Eshelby green's tensor (6x6 matrix)

% C0        : Stiffness tensor of the host material (6x6 matrix).
% alpha     : Alpha for a single inclusions.

I4 = eye(6);
S_0 = I4/C0; 
mu = C0(4,4)/2; 
kappa = C0(1,1) - (4/3)*mu; 
Poir = (3*kappa - 2*mu)/(2*(3*kappa + mu));

if alpha < 1
    q = (alpha/(1- alpha^(2))^(3/2))*(acos(alpha) - alpha*(1-alpha^(2))^(1/2));
    S_11 = (3/(8*(1-Poir)))*(alpha^(2)/(alpha^(2) - 1))+ (1/(4*(1-Poir)))*(1-2*Poir-9/(4*(alpha^(2)-1)))*q;
    S_33 = (1/(2*(1-Poir)))*(1-2*Poir + (3*alpha^(2)-1)/(alpha^(2)-1)-(1-2*Poir+3*alpha^(2)/(alpha^(2)-1))*q);
    S_12 = (1/(4*(1-Poir)))*(alpha^(2)/(2*(alpha^(2)-1))-(1-2*Poir + 3/(4*(alpha^(2)-1)))*q);
    S_13 = (1/(2*(1-Poir)))*(-(alpha^(2))/(alpha^(2)-1)+ 0.5*(3*alpha^(2)/(alpha^(2)-1)-(1-2*Poir))*q);

    S_31 = (1/(2*(1-Poir)))*(2*Poir - 1 - 1/(alpha^(2)-1)+ (1-2*Poir + 3/(2*(alpha^(2)-1)))*q);
    S_66 = (1/(4*(1-Poir)))*(alpha^(2)/(2*(alpha^(2)-1))+ (1-2*Poir -3/(4*(alpha^(2)-1)))*q);
    S_44 = (1/(4*(1-Poir)))*(1-2*Poir -(alpha^(2)+1)/(alpha^(2)-1)-0.5*(1-2*Poir-(3*(alpha^(2)+1))/(alpha^(2)-1))*q);
    S_22 = S_11;
    S_21 = S_12;
    S_23 = S_13;
    S_32 = S_31;
    S_55 = S_44;
    S_r = [S_11 S_12 S_13 0 0 0; S_21 S_22 S_23 0 0 0; S_31 S_32 S_33 0 0 0;
    0 0 0 2*S_44 0 0; 0 0 0 0 2*S_55 0; 0 0 0 0 0 2*S_66];
    G = -S_r*S_0;

elseif alpha == 1
    S_11 = (5*Poir - 1)/(15*(1 - Poir)) + (2*(4 - 5*Poir))/(15*(1 - Poir));
    S_12 = (5*Poir - 1)/(15*(1 - Poir));
    S_13 = (5*Poir - 1)/(15*(1 - Poir));
    S_31 = (5*Poir - 1)/(15*(1 - Poir));
    S_44 = (4 - 5*Poir)/(15*(1 - Poir));
    S_22 = S_11;
    S_33 = S_11;
    S_21 = S_12;
    S_32 = S_31;
    S_23 = S_13;
    S_55 = S_44;
    S_66 = S_44;
    S_r = [S_11 S_12 S_13 0 0 0; S_21 S_22 S_23 0 0 0; S_31 S_32 S_33 0 0 0;
            0 0 0 2*S_44 0 0; 0 0 0 0 2*S_55 0; 0 0 0 0 0 2*S_66];
    G = -S_r*S_0;
else
    G = 0;
end 
