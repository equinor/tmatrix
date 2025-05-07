function [ X_iso ] = iso_av_all( X, ctrl)

%Returns an multi dimesional matrix with isotropic elements. (6x6x(numbers of inclusions) matrix)  

% X : (6x6x(numbers of inclusions) matrix) if there is anisotropy in the
%      system, the first half of the matrices are the soon to become
%      isotropic. 
% ctrl : control parameter
% if ctrl = 0 then all the pore types are isotropic
% if ctrl(1) = 1 then all the pore types up to count(3)/2 are isotropic
% if ctrl(1) = 2 then all the pore types are anisotropic


% 09.03.2012
% Remy Agersborg 

count = size(X);
X_iso = zeros(6,6,count(3));
%if ctrl = 0 then all the pore types are isotropic
if ctrl(1) == 0 
    for j=1:(count(3))
%         X_temp(:,:)= X(:,:,j);
%         X_temp = iso_av(X_temp);
        X_iso(:,:,j) = iso_av(X(:,:,j));
    end;
%if ctrl(1) = 1 then all the pore types up to count(3)/2 are isotropic
elseif ctrl(1) == 1 
    for j=1:(count(3)/2)
%         X_temp(:,:)= X(:,:,j);
%         X_temp = iso_av(X_temp);
        X_iso(:,:,j) = iso_av(X(:,:,j));
    end;
    X_iso(:,:,(count(3)/2)+1:count(3)) = X(:,:,(count(3)/2)+1:count(3));
elseif ctrl(1) == 2
    X_iso(:,:,:) = X(:,:,:);
end    
end
