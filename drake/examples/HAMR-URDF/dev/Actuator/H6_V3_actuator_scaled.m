function actuator = H6_V3_actuator_scaled(t)

% Dimensions 
actuator.lext = 1;              % extension length in mm
actuator.l = 9.2;               % pzt length in mm
actuator.wn = 3.5;              % nominal width in mm
actuator.wr = 1.5;              % width ratio
actuator.tpzt = 135e-3;         % pzt-thickness (mm)   

%% PZT internal Stress as a function of field 

% For full model 
% actuator.f31_min = 14e-3;         %Pa*m/V (sigmoid min)
% actuator.f31_max = 29e-3;         %Pa*m/V (sigmoid max)
% actuator.field0 = 0.4e6;       %V/m (switching field)
% actuator.c = 1e-5;             %sigmoid width
% actuator.d = 69e-9;            %m/V (linear dependence of stress on field at high fields)

% for simple model 
actuator.f31_eff_b = 24e-3;         %Pa*mm/V (in blocked case)
actuator.f31_eff_f = 32e-3;         %Pa*mm/V (in free case)

%% PZT Young's modulus as a function of strain

% for full model 
actuator.Emin = 38.5e3;          %sigmoid min (Pa)
actuator.Emax = 81e3;            %sigmoid max (Pa)
actuator.eps0 = -0.00047;        %switching strain
actuator.a = 8e3;                  %sigmoid width

% for simple model 
actuator.Eave_eff = 40e3;        % young's modulus (Pa)

%% linear depedence of stress on strain 
% actuator.b = -230; 

%% CF propreties
actuator.tcf = t - 2*actuator.tpzt;       % carbon fiber thickness (m);
actuator.Ecf = 340e3;                     % young's modulus (Pa)