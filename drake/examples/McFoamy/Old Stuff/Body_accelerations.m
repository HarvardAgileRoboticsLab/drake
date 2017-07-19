function [qdd] = Body_accelerations(Vel,Fx,Fy,Fz,Mx,My,Mz,Quat,B_rate)

b = zeros(6,1);

%--------------------------------------------------------------------------
%                       Mass properties
%--------------------------------------------------------------------------
% Mass measured for flight on 18 April 2015
m = 0.484;
g = 9.81;

% Moments and products of inertia are from CAD model with config. v4
Ix  = 0.003922;
Iy  = 0.015940;
Iz  = 0.019340;
Ixz = 0.000441;
Ixy = 0.000303;
Iyz = -0.000030;

%--------------------------------------------------------------------------
%              Stability derivative - Added mass/inertia
%--------------------------------------------------------------------------

%-------------------------- X-force derivatives----------------------------

Xud = 0;   Xvd = 0;   Xwd = 0;   Xpd = 0;   Xqd = 0;   Xrd = 0;

%-------------------------- Y-force derivatives ---------------------------

Yud = 0;   Yvd = 0;   Ywd = 0;   Ypd = 0;   Yqd = 0;   Yrd = 0;

%-------------------------- Z-force derivatives ---------------------------

Zud = 0;   Zvd = 0;   Zwd = 0;   Zpd = 0;   Zqd = 0;   Zrd = 0;
% Zwd = -0.003912
%----------------------- Roll moment derivatives---------------------------

Lud = 0;   Lvd = 0;   Lwd = 0;   Lpd = 0;   Lqd = 0;   Lrd = 0;

%----------------------- Pitch Moment derivatives--------------------------

Mud = 0;   Mvd = 0;   Mwd = 0;   Mpd = 0;   Mqd = 0;   Mrd = 0;
% Mwd = -0.002073
%----------------------- Yaw Moment derivatives----------------------------

Nud = 0;   Nvd = 0;   Nwd = 0;   Npd = 0;   Nqd = 0;   Nrd = 0;

%--------------------------------------------------------------------------
%                       Defining Mass/Inertia matrix M
%--------------------------------------------------------------------------

M11 = m - Xud;
M12 = -Xvd;
M13 = -Xwd;
M14 = -Xpd;
M15 = -Xqd;
M16 = -Xrd;

M21 = -Yud;
M22 = m - Yvd;
M23 = -Ywd;
M24 = -Ypd;
M25 = -Yqd;
M26 = -Yrd;

M31 = -Zud;
M32 = -Zvd;
M33 = m - Zwd;
M34 = -Zpd;
M35 = -Zqd;
M36 = -Zrd;

M41 = -Lud;
M42 = -Lvd;
M43 = -Lwd;
M44 = Ix - Lpd;
M45 = -Ixy - Lqd; 
M46 = -Ixz - Lrd;

M51 = -Mud;
M52 = -Mvd;
M53 = -Mwd;
M54 = -Ixy - Mpd;
M55 = Iy - Mqd;
M56 = -Iyz - Mrd;

M61 = -Nud;
M62 = -Nvd;
M63 = -Nwd;
M64 = -Ixz - Npd;
M65 = -Iyz - Nqd;
M66 = Iz - Nrd;

M = [M11 M12 M13 M14 M15 M16; M21 M22 M23 M24 M25 M26; M31 M32 M33 M34 M35 M36;
 M41 M42 M43 M44 M45 M46; M51 M52 M53 M54 M55 M56; M61 M62 M63 M64 M65 M66];

e0 = Quat(1);
e1 = Quat(2);
e2 = Quat(3);
e3 = Quat(4);

% phi = Eul(1);
% the = Eul(2);
% psi = Eul(3);
% 
% cphi = cos(phi);
% cthe = cos(the);
% cpsi = cos(psi);
% sphi = sin(phi);
% sthe = sin(the);
% spsi = sin(psi);

u = Vel(1);
v = Vel(2);
w = Vel(3);
p = B_rate(1);
q = B_rate(2);
r = B_rate(3);

% To do : In the equations below, the masses and inertias should be
%         modified to include added mass/inertia - WHAT???
b(1) = Fx + m*g*2*(e1*e3 - e2*e0) + m*(r*v-q*w);
b(2) = Fy + m*g*2*(e2*e3 + e1*e0) + m*(p*w-r*u);
b(3) = Fz + m*g*(e3^2 + e0^2 - e1^2 - e2^2) + m*(q*u-p*v);
b(4) = Mx + (Iy-Iz)*(q*r) - Ixy*(p*r) - Iyz*(r^2-q^2)+ Ixz*(p*q);
b(5) = My + (Iz-Ix)*(p*r) + Ixy*(q*r) - Iyz*(p*q) - Ixz*(p^2-r^2);
b(6) = Mz + (Ix-Iy)*(p*q) - Ixy*(q^2-p^2) + Iyz*(p*r) - Ixz*(q*r);

qdd = M\b;