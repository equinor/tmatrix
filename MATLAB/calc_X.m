function [ X ] = calc_X(C0, td )
% Returns the X-tensor (6x6x(numbers of empty cavities) matrix) for more explanation see e.g.
% Agersborg et al. 2009 or "The effects of drained and undrained loading in
% visco-elsatic waves in rock-like composites" M. Jakobsen and T.A. % Johansen. 
%(2005). Int. J. Solids and Structures (42). p. 1597-1611 

% C0 :  Stiffness tensor of the host material (6x6 matrix)
% td :  Dry t-matrix tensors, (6x6x(numbers of empty cavities) matrix)

% 09.03.2012
% Remy Agersborg 

count = size(td);
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
X = zeros(6,6,count(3));

for j=1:count(3)
%     td_temp(:,:) = td(:,:,j);
%     X_temp = td_temp*S0*I2I2*S0*td_temp;
    X(:,:,j) =  td(:,:,j)*S0*I2I2*S0*td(:,:,j);%X_temp(:,:);
end;

