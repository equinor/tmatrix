function [Kd] = calc_Kd( C0, alpha)
% Returns the dry K-tensor for the inclusions specified by alpha
% Kd is a (6x6x(number of alphas)) matrix

% C0    : Stiffness tensor of the host material (6x6 matrix)
% alpha : vector of aspect ratios (1x (number of alphas) vector)

% Equations used can be found in:
% Agersborg (2007), phd thesis:
% https://bora.uib.no/handle/1956/2422

% 09.03.2012
% Remy Agersborg 

I4 = eye(6);

S0 = I4/C0;
L = length(alpha);

Kd=zeros(6,6,L);
% [a_un,pos] = unique(alpha,'first');
% 
% for nc = 1:length(a_un)
%     G  = Gtensor(C0, a_un(nc));
% %     Kd_temp = (I4 + G*C0)\S0;
%     Kd(:,:,[pos(nc) pos(nc)+L/2]) = [((I4 + G*C0)\S0);((I4 + G*C0)\S0)];
% %     Kd(:,:,pos(nc)+L/2) = ((I4 + G*C0)\S0);
%     
% end


for nc = 1:L
    G  = Gtensor(C0, alpha(nc));
    Kd(:,:,nc) =(I4 + G*C0)\S0;;
    
end

