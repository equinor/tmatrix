function out = T_matrix_porosity(mineral_property, fluid_property, phi_vector,scenario,  frequency, angle,per_inc_con,per_inc_ani)

%% T-matrix approach
%% by
%% Remy Agersborg
%% Bergen, 20.02.2012

%% The theory can be found in the papers and in the references therein:
%% Agersborg, R., Jakobsen, M., Ruud, B.O. and Johansen, T. A. 2007.
%% Effects of pore fluid pressure on the seismic response of a fractured carbonate reservoir.
%% Stud. Geophys. Geod., 51, 89-118.

%% Agersborg, R., Johansen, T. A. and Ruud, B.O. 2008.
%% Modelling reflection signatures of pore fluids and dual porosity in carbonate reservoirs.
%% Journal of Seismic Exploration, 17(1), 63-83.

%% Agersborg, R., Johansen, T. A., Jakobsen, M., Sothcott, J. and Best, A. 2008.
%% Effect of fluids and dual-pores systems on pressure-dependent velocities and attenuation in carbonates,
%% Geophysics, 73, No. 5, N35-N47.

%% Agersborg, R., Johansen, T. A., and Jakobsen, M. 2009.
%% Velocity variations in carbonate rocks due to dual porosity and wave-induced fluid flow.
%% Geophysical Prospecting, 57, 81-98.

%% All of the papers and a extended explanations of the involved equations
%% can be found in  Agersborg (2007), phd thesis: https://bora.uib.no/handle/1956/2422

% 09.03.2012
% Remy Agersborg

%% EDITED PART, ESTIMATE PARAMETERS IN MAIN FUNCITON
if scenario ==1 %Dual porosity, mostly rounded pores
    [alpha,v,tau]=deal(zeros(1,2));
    alpha(1) = 0.9;      v(1)= 0.9;   tau(1) = 10^-7;   % alpha: aspect ratio
    alpha(2) = 0.1;      v(2)= 0.1;    tau(2) = 10^-7;   % v : volume fraction
elseif scenario ==2 %Dual porosity, little rounded pores
    [alpha,v,tau]=deal(zeros(1,2));
    alpha(1) = 0.58;      v(1)= 0.85;    tau(1) = 10^-7;   % alpha: aspect ratio
    alpha(2) = 0.027;     v(2)= 0.15;    tau(2) = 10^-7;   % v : volume fraction
elseif scenario == 3 %Mixed pores
    [alpha,v,tau]=deal(zeros(1,3));
    alpha(1) = 0.9;      v(1)= 0.8;   tau(1) = 10^-7;   % alpha: aspect ratio
    alpha(2) = 0.1;      v(2)= 0.19;    tau(2) = 10^-7;   % v : volume fraction
    alpha(3) = 0.01;     v(3)= 0.01;    tau(3) = 10^-7;   % v : volume fraction
elseif scenario == 4 %Flat pores and cracks
    [alpha,v,tau]=deal(zeros(1,4));
    alpha(1) = 0.9;      v(1)= 0.689;   tau(1) = 10^-7;   % alpha: aspect ratio
    alpha(2) = 0.1;      v(2)= 0.3;    tau(2) = 10^-7;   % v : volume fraction
    alpha(3) = 0.01;     v(3)= 0.01;    tau(3) = 10^-7;   % v : volume fraction
    alpha(4) = 0.001;     v(4)= 0.001;    tau(4) = 10^-7;   % v : volume fraction
end

ctrl = zeros(1,2);
if (per_inc_ani ~= 0 && per_inc_ani ~= 1);
    ctrl(1) = 1;
end
if  per_inc_ani == 1;
    ctrl(1) = 2;
end

if per_inc_ani ~=0 && per_inc_ani ~=1 ;
    % divide the alphas, and v's into isotropy and anisotropy part (double the numbers of index where first part in the vector is isotropy)
    alpha = cat(2,alpha,alpha);
    v     = cat(2,(1-per_inc_ani)*v, per_inc_ani*v);
    tau   = cat(2,tau, tau);
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
%% Looping over porosity to calculate vectors of elastic properties
for ploop = 1:length(phi_vector)
    porosity = phi_vector(ploop);
    in_v.connected = initial_v.connected*porosity ;
    in_v.isolated = initial_v.isolated*porosity ;
    in_alpha = initial_alpha;
    
    in_v.connected(in_v.connected.*porosity>in_alpha.connected)=in_alpha.connected(in_v.connected.*porosity>in_alpha.connected)/2;
    in_v.isolated(in_v.isolated.*porosity>in_alpha.isolated)=in_alpha.isolated(in_v.isolated.*porosity>in_alpha.isolated)/2;
    
    %% ------- Fluid properties -------------------------------------------------
    % Permeability of the ROCK! MULTIPLICATE TO GET ON THE RIGHT SU UNITS
    K_r = fluid_property(ploop,3)*0.986923e-15;
    
    % Viscosity of the fluid
    eta_f = fluid_property(ploop,4)*1e-2; %MULTIPLICATE TO GET ON SU UNITS
    
    %Bulk modulus fluid
    if max(fluid_property(:,1))<1e9
        error('fluid_property should be a matrice with Bulk (Pa) Rho (kg/m^3) Perm (50) and eta')
    end
    kappa_f = fluid_property(ploop,1);
    rho_f = fluid_property(ploop,2);
    %% Creating reference matrix
    %Test that data are on correct format
    if max(mineral_property(:,3))<1000
        mineral_property(:,3) = mineral_property(:,3).*1e3;
    end
    if max(mineral_property(:,1:2)) < 1e9
        error('mineral properties must be a matrice with Bulk (Pa) Shear (Pa) and Rho (kg/m^3)')
    end
    rho_s = mineral_property(ploop,3);
    Vs_s = sqrt(mineral_property(ploop,2)./rho_s);
    Vp_s = sqrt((mineral_property(ploop,1)+4.*mineral_property(ploop,2)/3)./rho_s);
    
     c11 = rho_s.*Vp_s.^2;
    c44 = rho_s*Vs_s.^2;
    c12 = c11 - 2*c44;
    
    C0 = [
        c11 c12 c12 0     0     0;
        c12 c11 c12 0     0     0;
        c12 c12 c11 0     0     0;
        0   0   0   2*c44 0     0;
        0   0   0   0     2*c44 0;
        0   0   0   0     0     2*c44;
        ];
    
    
    %bulkmoddulus of hostmaterial of inclusions
    kappa = c11 - (4/3)*c44;
    %Green's function of interacting inclusions
    Gd = Gtensor(C0, 1);
    
    
    %Defining C1 tesor in the effective stiffness expression.
    C1 = zeros(6,6);
    if ctrl(2) ~= 2
        % Isolated part: calculated C1 tensor (sum over all the isolated
        % t-matrices and concentrations
        C1 = calc_isolated_part(C0,kappa_f,in_alpha.isolated, in_v.isolated, ctrl);
    end
    
    if ctrl(2) ~= 0
        % Connected part
        %Calculate dry properties:
        Kd = calc_Kd(C0,in_alpha.connected);
        td = calc_td(C0,Kd, in_alpha.connected);
        %iso averaging the isotropic porosity
        %NB: the anisotropic t-matrices are also included in td_bar
        td_bar = iso_av_all(td, ctrl);
        %Calculate fluid effect
        X = calc_X(C0,td);
        %iso averaging the isotropic porosity
        %NB: the anisotropic t-matrices are also included in td_bar
        X_bar = iso_av_all(X, ctrl);
        Kd_uuvv = calc_Kd_uuvv(Kd);
        gamma = 1 - kappa_f/kappa + kappa_f.*Kd_uuvv;
        %Frequency dependent Stiffness
        [C_eff ] = calc_C_eff_visco(Vs_s, K_r, eta_f, in_v.connected, gamma, tau, Kd_uuvv, kappa, kappa_f, C0, C1, td_bar, X_bar, Gd, frequency);
        
    end
    
    if ctrl(2) == 0
        I4 = eye(6);
        C_eff = C0 + C1/(I4 + Gd*C1);
    end
    
    
    
    %%%%%%%%%%%%%%
    
    % ctrl(2) = 0 :only isolated pores
    % ctrl(2) = 1 :both isolated and connected pores
    % ctrl(2) = 2 :only connected pores
    phi = 0;
    if ctrl(2) ~=  0
        phi = phi + sum(in_v.connected);
    end
    if ctrl(2) ~=  2
        phi = phi  + sum(in_v.isolated);
    end
    
    
    %Effective density
    rho_eff = phi*rho_f + (1-phi)*rho_s;
    if ctrl(1) == 0
        angle = 0;
    end
    
    [Vp, Vsv, Vsh ] = velocity_vti_angles( C_eff, rho_eff, angle );
    
    %output parameters
    out.Vp(ploop)  = Vp;
    out.Vsv(ploop)= Vsv;
    out.Vsh(ploop) = Vsh;
    out.frequency(ploop) = frequency;
    out.rho_eff(ploop) = rho_eff;
end

