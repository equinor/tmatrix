
function [ t ] = calc_t( td, Theta, X, Z, omega, gamma,tau,kappa_f)
% Returns the t-matrices (6x6x(numbers of connected pores)) of the
% connected pores 
% 
% td      : Dry t-matrix tensors, (6x6x(numbers of empty cavities) matrix).
% Theta   : Theta-tensor (6x6 matrix).
% X       : X-tensor (6x6x(numbers of empty cavities) matrix).
% Z       : Z-tensor (6x6x(numbers of empty cavities) matrix).
% omega   : frequency (2*pi*f).
% gamma   : gamma factor of all the inclusions (1x(numbers of empty
%           cavities) vector).
% tau     : Relaxation time constant (1x(numbers of empty cavities)
%           vector).
% kappa_f : Bulk modulus of the fluid.

% Equations and tensors used can be found in:
% Agersborg (2007), phd thesis:
% https://bora.uib.no/handle/1956/2422

% 09.03.2012
% Remy Agersborg 

a = ndims(td);
count = size(td);

if a == 2
   n = 1
else
    n = count(3);
end

t = zeros(6,6,n);

for j=1:n
    td_temp(:,:)= td(:,:,j);
    X_temp(:,:) = X(:,:,j);
    Z_temp(:,:) = Z(:,:,j);
    t(:,:,j) =  td(:,:,j) + (Theta*Z(:,:,j) + i*omega*tau(j)*kappa_f*X(:,:,j))/(1+i*omega*gamma(j)*tau(j));%td_temp + (Theta*Z_temp + i*omega*tau(j)*kappa_f*X_temp)/(1+i*omega*gamma(j)*tau(j));
end;


