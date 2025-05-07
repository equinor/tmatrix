function [ Kd_uuvv ] = calc_Kd_uuvv( Kd )
% Returns the sum of dry K_uuvv

% Kd : the dry K-tensor (6,6,(numbers of inclusions)) matrix

% Equations used can be found in:
% Agersborg (2007), phd thesis:
% https://bora.uib.no/handle/1956/2422

% 09.03.2012
% Remy Agersborg 

I2 = [1 1 1 0 0 0]';

a = ndims(Kd);
count = size(Kd);
if a == 2
    n = 1;
else
    n = count(3);
end
Kd_uuvv = zeros(1,n);

for j=1:n
%     Kd_temp(:,:) = Kd(:,:,j);
    Kd_uuvv(j) = I2'*Kd(:,:,j)*I2;
end;
