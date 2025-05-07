function [Vp, inv_Qp, Vs, inv_Qs] = calc_velocity_invQ(Ceff, rho_eff)

    dim_C = ndims(Ceff);
    count = size(Ceff);
    if dim_C == 2 
        count(3) = 1;
    end    
    
    Vp = zeros(count(3));
    inv_Qp = zeros(count(3));
    Vs = zeros(count(3));
    inv_Qs = zeros(count(3));
    
    for j = 1:count(3)
        %P-velocity & attenuation
        V2p    = Ceff(1,1,j)/rho_eff;               % velocity ^ 2
        sp     = real(1/sqrt(V2p));                 % real part of complex slowness
        Vp(j)  = 1/sp;                              % inverse of slowness (phase vel)
        %Qp     = abs(real(V2p)/imag(V2p));          % Quality factor
        %inv_Qp(j) = 1/Qp;
     
    
        %S-velocity & attenuation
        V2s    = Ceff(4,4,j)/(2*rho_eff);
        ss     = real(1/sqrt(V2s));
        Vs(j)     = 1/ss;
        %Qs     = abs(real(V2s)/imag(V2s));
        %inv_Qs(j) = 1/Qs;
    end
    
    %%%%inv_Qp _Qs set to zero
    inv_Qp =0;
    inv_Qs =0;
end
