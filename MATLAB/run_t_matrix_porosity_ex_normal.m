%% Defining parameters. Initializing the model ---->
loglength = 1e3;
% r1 = 0.8 + (1.15 - 0.8).*rand(loglength,1); %ones(loglength,1);
% r2 = 0.8 + (1.15 - 0.8).*rand(loglength,1); %ones(loglength,1);
r1= ones(loglength,1);
r2 = ones(loglength,1);
phi_vector = (linspace(0.001,0.33,loglength))';%(rand(loglength,1)./2;
%------- Mineral properties -----------------------------------------------
Vp = 6260;                                    % P-velocity (m/s)
Vs = 3240;                                    % S-velocity (m/s)
rho= 2710;                                    % Density (kg/m3)

Gs =Vs.^2.*rho;                               % Use moduli instead of V
Ks = rho.*Vp.^2 - 4*Gs/3;
mineral_property = [Ks.*r1, Gs.*r2, rho.*(r1+r2)/2];

%------- Fluid properties -------------------------------------------------

rho_f = 1000;                                   % Fluid density (kg/m3)
V_f = 1500;                                     % Fluid velocity (m/s)

K_r = 50;                                       % Permeability (mDa) OF ROCK !! NOT FLUID RELATED.
eta_f = 1;                                      % Viscocity (cP)

kappa_f = rho_f.*V_f.^2;                            % Use moduli instead of V

fluid_property = [kappa_f.*r1, rho_f.*r2, K_r.*r1, eta_f.*r1];

%----------------------------- Anisotropy  -----------------------
per_inc_ani = 0.4;                   %percentage (0-1) of inclusions which is anisotropic
angle_of_sym_plane = 45;            %Angle (degree) of symmetry plane =0 if HTI media and =90 if VTI media

%-------------------- Connectivity of pores -----------------------
per_inc_con = 0.5 ;                   %percentage (0-1) of the inclusions which is connected

frequency = 50;
scenario = 3;
tic
out = T_matrix_porosity(mineral_property, fluid_property, phi_vector, scenario,  frequency, angle_of_sym_plane,per_inc_con,per_inc_ani);
toc
Vp  = out.Vp;
Vsh = out.Vsh;
Vsv = out.Vsv;
rho_eff = out.rho_eff;

phi_tot = phi_vector;

a=real(Vsv);
b = find(Vsv==a);
ai = rho_eff(b).*Vp(b)./1e6;
vpvs = Vp(b)./Vsv(b);

figure(1);hold on
%subplot(3,1,1);
plot(phi_tot, Vp,'m.','markersize',15);
hold on;
ylabel('Vp (m/s)');xlabel('Porosity');
hold on;

figure(2);hold on
hold on
%subplot(3,1,2);
plot(phi_tot, Vsv,'m.','markersize',15);
hold on;
ylabel('Vsv (m/s)');xlabel('Porosity');

figure(3);hold on
% plot(phi_tot, Vsh./Vsv,'r.','markersize',15);
plot(phi_tot, Vsh,'m.','markersize',15);
ylabel('Vsh (m/s)');xlabel('Porosity');
hold on;

rho = por2dens(phi_tot,0);
AI = rho'.*Vp./1e6;
VpVs = Vp./Vsv;

indx = find(AI>0 & VpVs>0);

