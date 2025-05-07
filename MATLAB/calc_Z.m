function [ Z ] = calc_Z( C0,td, omega, gamma, v, tau )

% Returns the Z tensor (6x6x(numbers of empty cavities) matrix) for more explanation see e.g.
% Agersborg et al. 2009 or "The effects of drained and undrained loading in
% visco-elsatic waves in rock-like composites" M. Jakobsen and T.A. % Johansen. 
%(2005). Int. J. Solids and Structures (42). p. 1597-1611 

% C0    : Stiffness tensor of the host material (6x6 matrix)
% td    : Dry t-matrix tensors, (6x6x(numbers of empty cavities) matrix)
% omega : frequency (2*pi*f)
% gamma : gamma factor of all the inclusions (1x(numbers of empty cavities)
%         vector) 
% v     :  concentration of all the empty cavities (1x(numbers of empty cavities)
%          vector) 
% tau   :  Relaxation time constant (1x(numbers of empty cavities)
%          vector)

% 09.03.2012
% Remy Agersborg 

I2I2 = [
1 1 1 0 0 0
1 1 1 0 0 0
1 1 1 0 0 0
0 0 0 0 0 0
0 0 0 0 0 0
0 0 0 0 0 0
];

I4 = eye(6);
S0 = I4/C0;

sum=0;
Z = zeros(6,6,length(v));
for j=1:length(v)
%     td_temp(:,:)=td(:,:,j);
    sum = sum + v(j)*td(:,:,j)/(1+i*omega*gamma(j)*tau(j));
end;
for j=1:length(v)
%     td_temp(:,:)=td(:,:,j);
%     Z_temp= td_temp*S0*I2I2*S0*sum;
    Z(:,:,j) = td(:,:,j)*S0*I2I2*S0*sum; %Z_temp(:,:);
end;
