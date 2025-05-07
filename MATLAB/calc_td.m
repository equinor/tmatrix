function [td] = calc_td( C0, Kd, alpha)
% Returns the dry t-matrix tensors (6x6x(numbers of empty cavities) matrix)

%C0    : Stiffness tensor of the host material (6x6 matrix)
%Kd    : Dry K tensor of all the empty cavities (6x6x(numbers of empty
%        cavities) matrix) see Agersborg et al. 2009 for explanation. 
%alpha : aspect ratios of all the empty cavities (1x(numbers of empty
%        cavities) vector)

% 09.03.2012
% Remy Agersborg 

I4 = eye(6);
L = length(alpha);
td = zeros(6,6,L);

for nc = 1:L
    G  = Gtensor(C0, alpha(nc));
%     Kd_temp = Kd(:,:,nc);
%     td_temp = G\(Kd_temp*C0 - I4);
    td(:,:,nc) =  G\(Kd(:,:,nc)*C0 - I4);% td_temp(:,:);
end

