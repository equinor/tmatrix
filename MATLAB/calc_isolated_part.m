function [ C1 ] = calc_isolated_part( C0, kappa_f, alpha, v, ctrl )
%  Returns the first order correction tensor: sum of the concentrations and
%  t-matrices of the isolated porosity (6x6 matrix).

% C0      : Stiffness tensor of the host material
% kappa_f : bulk modulus of the fluid
% alpha   : aspect ratios of all the inclusions (1xnumber of inclusions) vector)
% v       : concentration of all the inclusions (1xnumber of inclusions) vector)
% ctrl    : control vector 
% ctrl(1) = 0 : 100% isotropic
% ctrl(1) = 1 : mixed isotropic and anisotropic porosity
% ctrl(1) = 2 : 100% anisotropic porosity

% Equations used can be found in:
% Agersborg (2007), phd thesis:
% https://bora.uib.no/handle/1956/2422

% 09.03.2012
% Remy Agersborg 

%set C^n (Stiffness of the material inside the inclusion)
Cn = zeros(6,6);
Cn(1:3,1:3) = kappa_f;

I4 = eye(6);
C1 = zeros(6,6);
if ctrl(1) ~= 1
    %if there is only isotropic or anisortopic inclusions
    for j=1:length(alpha)
        G  = Gtensor(C0, alpha(j));
        t =  (Cn - C0)/(I4 - G*(Cn -C0));
        if ctrl(1) ~= 2  
            t = iso_av(t);
        end;
        C1 = C1 + v(j)*t;
    end
else
    %Isotropic part
    for j=1:length(alpha)/2
        G  = Gtensor(C0, alpha(j));
        t =  (Cn - C0)/(I4 - G*(Cn -C0));
        t = iso_av(t);
        C1 = C1 + v(j)*t;
    end
    %anisotropic part
    for j=length(alpha)/2+1:length(alpha)
        G  = Gtensor(C0, alpha(j));
        t =  (Cn - C0)/(I4 - G*(Cn -C0));
        C1 = C1 + v(j)*t;
    end
end

end
        

