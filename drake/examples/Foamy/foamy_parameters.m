function p = foamy_parameters()

% ----- Model Parameters for McFoamy RC Airplane ----- %

g = 9.81; %Gravitational acceleration (m/s^2)
rho = 1.2; %Air density at 20C (kg/m^3)
m = .484; %Mass of plane (kg)

% Inertias
Jx  = 0.003922;
Jy  = 0.015940;
Jz  = 0.019340;
Jxz = 0.000441;
Jxy = 0.000303;
Jyz = -0.000030;
J = [Jx  Jxy Jxz;
     Jxy Jy  Jyz;
     Jxz Jyz Jz];
Jinv = inv(J);

% Propeller
g_thr = .02; %maps control input to Newtons of thrust (approximate)
r_prop = 12.7/100;
D_prop = r_prop*2;
A_prop = pi*r_prop*r_prop;
Jprop = .015*r_prop*r_prop/3; %motor + prop inertia (kg*m^2) %TODO: Measure this
%Coefficients from Selig for APC Slow Flyer 10x4.7
%These are realy functions of advance ratio, but approximating as constants for now
Ct = 0.1; %Thrust coefficient
Cp = 0.05; %Power coefficient
Cq = Cp/(2*pi); %Torque coefficient

% All lifting surfaces are modeled as unsweapt tapered wings
b = 86.4/100; %wing span (m)
l_in = 20/100; %inboard wing length covered by propwash (m)
cr = 26/100; %root chord (m)
ct = 15.2/100; %tip chord (m)
cm = (ct + cr)/2; %mean wing chord (m)
S = b*cm; %planform area of wing (m^2)
S_in = 2*l_in*cr;
S_out = S-S_in;
Ra = b^2/S; %wing aspect ratio (dimensionless)
Rt = ct/cr; %wing taper ratio (dimensionless)
r_ail = (b/6)*(1+2*Rt)/(1+Rt); %aileron moment arm (m)

ep_ail = 0.8; %flap effectiveness (Phillips P.41)
trim_ail = 127; %control input for zero deflection
g_ail = (45*pi/180)/255; %maps control input to deflection angle %TODO: Calibrate

b_elev = 18.2/100; %elevator span (m)
cr_elev = 15.2/100; %elevator root chord (m)
ct_elev = 14/100; %elevator tip chord (m)
cm_elev = (ct_elev + cr_elev)/2; %mean elevator chord (m)
S_elev = b_elev*cm_elev; %planform area of elevator (m^2)
Ra_elev = b_elev^2/S_elev; %elevator aspect ratio (dimensionless)
r_elev = 45/100; %elevator moment arm (m)

ep_elev = 0.9; %flap effectiveness (Phillips P.41)
trim_elev = 127; %control input for zero deflection
g_elev = (45*pi/180)/255; %maps control input to deflection angle %TODO: Calibrate

b_rud = 21.6/100; %rudder span (m)
cr_rud = 20.4/100; %rudder root chord (m)
ct_rud = 12.9/100; %rudder tip chord (m)
cm_rud = (ct_rud + cr_rud)/2; %mean rudder chord (m)
S_rud = b_rud*cm_rud; %planform area of rudder (m^2)
Ra_rud = b_rud^2/S_rud; %rudder aspect ratio (dimensionless)
r_rud = 48/100; %rudder moment arm (m)
z_rud = 3/100; %height of rudder center of pressure (m)

ep_rud = 0.9; %flap effectiveness (Phillips P.41)
trim_rud = 127; %control input for zero deflection
g_rud = (45*pi/180)/255; %maps from control input to deflection angle %TODO: Calibrate

%Lift curve polynomial fit
Clcoef = [-9.781885297556400 38.779513049043175 -52.388499489940138 19.266141214863080 15.435976905745736 -13.127972418509980 -1.155316115022734 3.634063117174400 -0.000000000000001]';
Cdcoef = [-0.353585618276247 3.607550808703421 -10.272069825351215 -4.489225907857385 -2.746985301074068 3.480420330498847 0.085904634206004 0.063691497636087]';

% --- Pack everything into a struct --- %
p = struct();
p.g = g;
p.rho = rho;
p.m = m;
p.Jx = Jx;
p.Jy = Jy;
p.Jz = Jz;
p.J = J;
p.Jinv = Jinv;
p.g_thr = g_thr;
p.r_prop = r_prop;
p.D_prop = D_prop;
p.A_prop = A_prop;
p.Jprop = Jprop;
p.Ct = Ct;
p.Cp = Cp;
p.Cq = Cq;
p.b = b;
p.l_in = l_in;
p.cr = cr;
p.ct = ct;
p.cm = cm;
p.S = S;
p.S_in = S_in;
p.S_out = S_out;
p.Ra = Ra;
p.Rt = Rt;
p.r_ail = r_ail;
p.ep_ail = ep_ail;
p.trim_ail = trim_ail;
p.g_ail = g_ail;
p.b_elev = b_elev;
p.cr_elev = cr_elev;
p.cm_elev = cm_elev;
p.S_elev = S_elev;
p.Ra_elev = Ra_elev;
p.r_elev = r_elev;
p.ep_elev = ep_elev;
p.trim_elev = trim_elev;
p.g_elev = g_elev;
p.b_rud = b_rud;
p.cr_rud = cr_rud;
p.ct_rud = ct_rud;
p.cm_rud = cm_rud;
p.S_rud = S_rud;
p.Ra_rud = Ra_rud;
p.r_rud = r_rud;
p.z_rud = z_rud;
p.ep_rud = ep_rud;
p.trim_rud = trim_rud;
p.g_rud = g_rud;
p.Clcoef = Clcoef;
p.Cdcoef = Cdcoef;

end
