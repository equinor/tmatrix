function [ Vp_out, Vsv_out, Vsh_out ] = velocity_vti_angles( Ceff, rho_eff, angle )
%Returns the P-velocity and  S-velocities.
%
%Ceff   :	Effective stiffness tensor (6x6 matrix)
%rho_eff:	Effective density
%angle	:   The angle between the wave vector and the axis of symmetry.


%Vp90
V2p= (Ceff(1,1))/rho_eff;
sp= real(1/(sqrt(V2p)));
Vp90 = 1/sp;

%Vp0
V2p= (Ceff(3,3))/rho_eff;
sp= real(1/(sqrt(V2p)));
Vp0 = 1/sp;

%Vp45
m45 = ((Ceff(1,1) - (Ceff(4,4))/2)*0.5 - (Ceff(3,3)-(Ceff(4,4))/2)*0.5 )^2  + (Ceff(1,1) - (Ceff(4,4))/2)^2;
m45 = sqrt(m45);
V2p = (Ceff(1,1)*0.5 + Ceff(3,3)*0.5 +  Ceff(4,4)/2 + m45)/(2*rho_eff);
sp= real(1/(sqrt(V2p)));
Vp45 = 1/sp;

%Vsh90
V2s= Ceff(6,6)/(2*rho_eff);
ss= real(1/(sqrt(V2s)));
Vs90 = 1/ss;

%Vsh0
V2s= Ceff(4,4)/(2.d0*rho_eff);
ss= real(1.d0/(sqrt(V2s)));
Vs0 = 1/ss          ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C11 = rho_eff*Vp90^2;
C33 = rho_eff*Vp0^2;
C44 = rho_eff*Vs0^2;
C66 = rho_eff*Vs90^2;
C13 = -C44 + sqrt( 4*(rho_eff^2)*(Vp45^4) - 2*rho_eff*(Vp45^2)*(C11 + C33 + 2*C44)+ (C11 + C44)*(C33 + C44));

rad_angle= (angle*pi)/180;

m_real = ((C11 -C44)*(sin(rad_angle)^2) - (C33-C44)*(cos(rad_angle)^2))^2 + ((C13 + C44)^2)*(sin(2*rad_angle)^2);

Vp_out  = sqrt(C11*(sin(rad_angle)^2) + C33*(cos(rad_angle)^2) + C44 + sqrt(m_real))*sqrt(1/(2*rho_eff));

Vsv_out = sqrt(C11*(sin(rad_angle)^2)   + C33*(cos(rad_angle)^2) + C44 - sqrt(m_real))*sqrt(1/(2*rho_eff));

Vsh_out = sqrt((C66*(sin(rad_angle)^2) + C44*(cos(rad_angle)^2))/rho_eff);

end