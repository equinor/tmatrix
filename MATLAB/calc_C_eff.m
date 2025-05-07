function [C_eff] = calc_C_eff (C0, C1, Gd, t, v)
% Find t-matrix for each inclusion
% Equation  4 (page 222) Jakobsen et al. 2003 (The acoustic signature of fluid flow in complex
% porous media)

%Returns the effective stiffness tensor C* (6x6 matrix) calculated from the t-matrices t(r) (Eq. 1). 
%C0	: Stiffness tensor of the host of the inclusion (6x6 matrix)
%C1	: Sum of the concentration and t-matrices (6x6 matrix)
%Gd	: Correlation function (6x6 matrix)
%t	: t-matrices of the different inclusions (6x6xr matrix)
%v	: Concentration of the inclusions (1xr vector)

% r : number of inclusions

% Equation can be found in:
% Agersborg (2007), phd thesis:
% https://bora.uib.no/handle/1956/2422


% 09.03.2012
% Remy Agersborg 


I4 = eye(6);
count = size(t);
for j = 1:count(3)
%     Cj = v(j)*t(:,:,j);
    C1 = C1 + v(j)*t(:,:,j);
end

C_eff = C0 + C1/(I4 + Gd*C1);


end