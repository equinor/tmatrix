
%% Defining parameters. Initializing the model ---->

%------- Mineral properties -----------------------------------------------
Vp = 6260;                                    % P-velocity (m/s)
Vs = 3240;                                    % S-velocity (m/s)
rho= 2710;                                    % Density (kg/m3)

mineral_property = [Vp, Vs, rho];
%------- Fluid properties -------------------------------------------------

rho_f = 1000;                                   % Fluid density (kg/m3)
V_f = 1500;                                     % Fluid velocity (m/s)

K_r = 50;                                       % Permeability (mDa)
eta_f = 1;                                      % Viscocity (cP)

fluid_property = [V_f, rho_f, K_r, eta_f];


%--------- Porosity ---------------------------------------------------------

porosity_range = [0.01, 0.45];    %range of porosity [initial porosity, final porosity]
dphi           = 0.01;             %steps of increasing the porosity

%----------------------------- Anisotropy  -----------------------
per_inc_ani = 0.4;                   %percentage (0-1) of inclusions which is anisotropic
angle_of_sym_plane = 45;            %Angle (degree) of symmetry plane =0 if HTI media and =90 if VTI media

%-------------------- Connectivity of pores -----------------------
per_inc_con = 0.5 ;                   %percentage (0-1) of the inclusions which is connected
%with respect to fluid flow.
%NB: Not permeabillity



%------- Porosity in the host to reduce mineral properties ----------------
% Soft porosity example: max total porosity=0.45
alpha(1) = 0.6;      v(1)= 0.899;    tau(1) = 10^-7;   % alpha: aspect ratio
alpha(2) = 0.1;      v(2)= 0.1;      tau(2) = 10^-7;   % v : volume fraction
alpha(3) = 0.001;    v(3)= 0.001;    tau(3) = 10^-7;

% Check if v(j) > alpha(j)for maximum porosity. If true, set v(j) = alpha(j)/2 to make sure
% the numbers of inclusions in the system is not violating the
% approximations for effective medium theories.
v = v*porosity_range(2);
[alpha, v] = check_alpha_concentration( alpha, v );
v = v./porosity_range(2);
%%Sort the porosity according to parameters above

% control parameter:
% ctrl(1) = 0 : 100% isotropy
% ctrl(1) = 1 : mixed isotropy and anisotropy
% ctrl(1) = 2 : 100% anisotropy

ctrl = zeros(1,2);
if (per_inc_ani ~= 0 && per_inc_ani ~= 1);
    ctrl(1) = 1;
end
if  per_inc_ani == 1;
    ctrl(1) = 2;
end

if per_inc_ani ~=0 && per_inc_ani ~=1 ;
    % divide the alphas, and v's into isotropy and anisotropy part (double the numbers of index where first part in the vector is isotropy)
    [alpha, v, tau] = sort_pore_types(alpha, v, tau, per_inc_ani);
end

%divide the porosity into isolated and connected part
if per_inc_con ~= 0  && per_inc_con ~= 1
    ctrl(2) = 1;
end
if per_inc_con == 1
    ctrl(2) =2;
end

%Finalizing the initial inclusions (porosity)
initial_alpha.connected = alpha;
initial_alpha.isolated  = alpha;
initial_v.connected = v.*per_inc_con;
initial_v.isolated = v.*(1-per_inc_con);
%% end sorting the porosity


%Set the porosity
if porosity_range(1) == porosity_range(2)
    phi_vector = porosity_range(1);
else
    phi_vector = porosity_range(1):dphi:porosity_range(2);
end


%%%%% End - this should never be inside a loop

frequency = 50;

%define Vp, Vsv, Vsh for 100 different porosities
count_porosity = size(phi_vector);

Vp = zeros(length(phi_vector));
Vsh = zeros(length(phi_vector));
Vsv = zeros(length(phi_vector));
phi_tot = zeros(length(phi_vector));

for j = 1:(length(phi_vector))
    in_v.connected = initial_v.connected*(phi_vector(j)) ;
    in_v.isolated = initial_v.isolated*(phi_vector(j)) ;
    in_alpha = initial_alpha;
    
    %     out = T_matrix_porosity(mineral_property, fluid_property, in_alpha, in_v, tau,  frequency, angle_of_sym_plane,ctrl);
    out = T_matrix_porosity(mineral_property, fluid_property, phi_vector(j), alpha, v, tau,  frequency, angle_of_sym_plane,per_inc_con,per_inc_ani);
    
    Vp(j,:)  = out.Vp;
    Vsh(j,:) = out.Vsh;
    Vsv(j,:) = out.Vsv;
    phi_tot(j,:) = out.phi_tot;
    
end

%This should be removed. Only for the example.
%Plotting the velocities
figure(1)
plot(phi_tot, Vp,'-k','linewidth',3);
hold on;
ylabel('Vp (m/s)');xlabel('Porosity');
hold on;

figure(2)
hold on;
%subplot(3,1,2);
plot(phi_tot, Vsv,'-k','linewidth',3);
plot(phi_tot, Vsh,'-r','linewidth',3);
ylabel('Vsv and Vsh (m/s)');xlabel('Porosity');

figure(3)
plot(phi_tot, Vsh./Vsv,'-m','linewidth',3);
ylabel('Vsh/Vsv');xlabel('Porosity');
hold on;


