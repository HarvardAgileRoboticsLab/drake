function [Fx,Fy,Fz,Mx,My,Mz] = Aerodynamics(Vel,CtrlDef,ThrForce,ThrMom,GyroMom,VProp_Axial,VProp_Swirl,V_wind,Geometry,CG,B_rate)
rho = 1.225;    % Density of air
nu = 1.56e-05;  % Kinematic viscosity of air
Re_crit = 3.e5;

%--------------------------------------------------------------------------
%                Geometric properties of components
%--------------------------------------------------------------------------
% AIRCRAFT CG
% Measured w.r.t the nose (from propeller plane)
xCG = CG(1)*1e-3;
yCG = CG(2)*1e-3;
zCG = CG(3)*1e-3;

% THRUSTER
% Thrust is assumed to act in x direction only
xp = 0 - xCG;
yp = 0 - yCG;
zp = 0 - zCG;

% WING GEOMETRY
Nw = 7; %Geometry(1,1);
Cw0 = Geometry(2,1)*1e-3;
Cwe = Geometry(3,1)*1e-3;
Bw = Geometry(4,1)*1e-3;
AR_w = (Bw)/(0.5*(Cw0 + Cwe));  % Aspect ratio

bw = Geometry(1,5:Nw+4)*1e-3;
Aw = Geometry(2,5:Nw+4)*1e-6;
Cw = Geometry(3,5:Nw+4)*1e-3;
Cwf = Geometry(4,5:Nw+4)*1e-3;
xw = Geometry(5,5:Nw+4)*1e-3- xCG;
yws = Geometry(6,5:Nw+4)*1e-3 - yCG;
ywp = -yws;
zw = Geometry(7,5:Nw+4)*1e-3 - zCG;

% TAIL GEOMETRY
Nt = 3; %Geometry(1,2);
Ct0 = Geometry(2,2)*1e-3;
Cte = Geometry(3,2)*1e-3;
Bt = Geometry(4,2)*1e-3;
AR_t = (Bt)/(0.5*(Ct0 + Cte));  % Aspect ratio

bt = Geometry(1,Nw+5:Nw+Nt+4)*1e-3;
At = Geometry(2,Nw+5:Nw+Nt+4)*1e-6;
Ct = Geometry(3,Nw+5:Nw+Nt+4)*1e-3;
Ctf = Geometry(4,Nw+5:Nw+Nt+4)*1e-3;
xt = Geometry(5,Nw+5:Nw+Nt+4)*1e-3 - xCG;
yts = Geometry(6,Nw+5:Nw+Nt+4)*1e-3 - yCG;
ytp = -yts;
zt = Geometry(7,Nw+5:Nw+Nt+4)*1e-3 - zCG;
 
% RUDDER GEOMETRY
Nr = 4; %Geometry(1,3);
Cr0 = Geometry(2,3)*1e-3;
Cre = Geometry(3,3)*1e-3;
Br = Geometry(4,3)*1e-3;
AR_r = Br/(0.5*(Cr0 + Cre));  % Aspect ratio

br = Geometry(1,Nw+Nt+5:Nw+Nt+Nr+4)*1e-3;
Ar = Geometry(2,Nw+Nt+5:Nw+Nt+Nr+4)*1e-6;
Cr = Geometry(3,Nw+Nt+5:Nw+Nt+Nr+4)*1e-3;
Crf = Geometry(4,Nw+Nt+5:Nw+Nt+Nr+4)*1e-3;
xr = Geometry(5,Nw+Nt+5:Nw+Nt+Nr+4)*1e-3 - xCG;
yr = Geometry(6,Nw+Nt+5:Nw+Nt+Nr+4)*1e-3 - yCG;
zr = Geometry(7,Nw+Nt+5:Nw+Nt+Nr+4)*1e-3 - zCG;

% BODY GEOMETRY
NB = 4; %Geometry(1,4);
CB0 = Geometry(2,4)*1e-3;
CBe = Geometry(3,4)*1e-3;
BB = Geometry(4,4)*1e-3;
AR_B = BB/(0.5*(CB0 + CBe));  % Aspect ratio

bB = Geometry(1,Nw+Nt+Nr+5:Nw+Nt+Nr+NB+4)*1e-3;
AB = Geometry(2,Nw+Nt+Nr+5:Nw+Nt+Nr+NB+4)*1e-6;
CB = Geometry(3,Nw+Nt+Nr+5:Nw+Nt+Nr+NB+4)*1e-3;
CBf = Geometry(4,Nw+Nt+Nr+5:Nw+Nt+Nr+NB+4)*1e-3;
xB = Geometry(5,Nw+Nt+Nr+5:Nw+Nt+Nr+NB+4)*1e-3 - xCG;
yB = Geometry(6,Nw+Nt+Nr+5:Nw+Nt+Nr+NB+4)*1e-3 - yCG;
zB = Geometry(7,Nw+Nt+Nr+5:Nw+Nt+Nr+NB+4)*1e-3 - zCG;

% WING TAIL GEOMETRY
GTermTS_ws = [0.016084,0.023292,0.02174,0.024746,0.019789,0.008994,0.009452;...
    0.015031,0.022185,0.021276,0.024887,0.020451,0.009449,0.010097;...
    0.013886,0.02128,0.021401,0.026499,0.022854,0.010896,0.012031;...
    0.012501,0.019888,0.021006,0.027634,0.025349,0.012542,0.014333];

GTermTS_wp = [0.016028,0.022902,0.021162,0.023799,0.01885,0.008538,0.008974;...
    0.014821,0.020897,0.019042,0.021261,0.016865,0.007657,0.007941;...
    0.013364,0.018283,0.01618,0.017651,0.013859,0.006286,0.006349;...
    0.0118,0.015694,0.013562,0.014531,0.011319,0.005122,0.00507];

GTermTP_ws = GTermTS_wp;
GTermTP_wp = GTermTS_ws;

delGTermTS_ws = [9.679512,-2.098619,-0.735873,-0.430774,-0.263924,-0.179803,-0.171809,-0.149025,-0.141583;...
    2.063419,-7.128762,-0.995697,-0.519464,-0.302859,-0.202192,-0.194671,-0.166454,-0.159319;...
    0.821159,1.609437,-4.291581,-0.908826,-0.427058,-0.264231,-0.257272,-0.212266,-0.205356;...
    0.470922,0.674116,1.671211,-3.911107,-0.709318,-0.371399,-0.365428,-0.28553,-0.279139];

delGTermTS_wp = [9.679512,1.457408,0.631687,0.388822,0.244562,0.168827,0.160982,0.140645,0.133319;...
    2.276448,0.970612,0.51223,0.335857,0.219725,0.15519,0.147915,0.130546,0.123809;...
    0.821159,0.542403,0.351353,0.251192,0.174517,0.127867,0.121256,0.109256,0.10323;...
    0.470922,0.355818,0.254443,0.191957,0.139323,0.105163,0.099164,0.091074,0.085642];

delGTermTP_ws = - delGTermTS_wp;
delGTermTP_wp = - delGTermTS_ws;

% WING RUDDER GEOMETRY
GTermR_w = [0.01091,0.015735,0.014915,0.017615,0.01486,0.00711,0.00774;...
    0.011013,0.015935,0.015057,0.017762,0.01494,0.007135,0.00779;...
    0.010909,0.015759,0.014912,0.017611,0.014856,0.007088,0.007754;...
    0.010155,0.014703,0.013956,0.016596,0.014143,0.006827,0.007436];

ThetaR_w = [3.22,3.24,3.26,3.28,3.3,3.33,3.64;...
    -2.84,-2.86,-2.87,-2.89,-2.92,-2.93,-3.21;...
    -8.57,-8.61,-8.66,-8.72,-8.79,-8.84,-9.67;...
    -13.53,-13.59,-13.67,-13.76,-13.86,-13.94,-15.19]*pi/180;

%--------------------------------------------------------------------------
%                       Aerodynamic Constants
%--------------------------------------------------------------------------
Cd0 = 0.02;
Cd90 = 1.98;
alp0 = 0;
% alpStallP = 10*pi/180;
% alpStallN = -10*pi/180;
% HighAlpStart = 22*pi/180;
% HighAlpEnd = 158*pi/180;
% StallMdl = 'FullStall';

%--------------------------------------------------------------------------
%                       Wind disturbances
%--------------------------------------------------------------------------
% Wind disturbance in body frame
Vw_u = V_wind(1);
Vw_v = V_wind(2);
Vw_w = V_wind(3);
% Vw_p = V_wind(4);
% Vw_q = V_wind(5);
% Vw_r = V_wind(6);
Vw_p = 0;
Vw_q = 0;
Vw_r = 0;

%--------------------------------------------------------------------------
%     Body axis components of vehicle velocity relative to the air
%--------------------------------------------------------------------------
% Assuming wind velocity is in body frame
u = Vel(1)- Vw_u;%waqas has +
v = Vel(2)- Vw_v;%waqas has +
w = Vel(3)- Vw_w;%waqas has +

p = B_rate(1)+ Vw_p;
q = B_rate(2)+ Vw_q;
r = B_rate(3)+ Vw_r;

%--------------------------------------------------------------------------
%              Slipstream velocities at reference points
%--------------------------------------------------------------------------
Vp_w = VProp_Axial(1:Nw);
Vp_t = VProp_Axial(Nw+1:Nw+Nt);
Vp_r = VProp_Axial(Nw+Nt+1:Nw+Nt+Nr);
Vp_B = VProp_Axial(Nw+Nt+Nr+1:Nw+Nt+Nr+NB);

Vp_w_swirl = VProp_Swirl(1:Nw);
Vp_t_swirl = VProp_Swirl(Nw+1:Nw+Nt);
Vp_r_swirl = VProp_Swirl(Nw+Nt+1:Nw+Nt+Nr);
Vp_B_swirl = VProp_Swirl(Nw+Nt+Nr+1:Nw+Nt+Nr+NB);

%--------------------------------------------------------------------------
%                             Control Inputs 
%--------------------------------------------------------------------------
RAilDef = CtrlDef(1);
LAilDef = CtrlDef(2);
ElevDef = CtrlDef(3);
RudDef = CtrlDef(4);

%--------------------------------------------------------------------------
%                   Wing, Tail, Rudder and Body velocities
%--------------------------------------------------------------------------
vws_xcomp = u + q*zw - r*yws + Vp_w;
vws_ycomp = v + r*xw - p*zw;
vws_zcomp = w + p*yws - q*xw - Vp_w_swirl;
vws = sqrt(vws_xcomp.^2 + vws_ycomp.^2 + vws_zcomp.^2);
vws_xz = sqrt(vws_xcomp.^2 + vws_zcomp.^2);

vwp_xcomp = u + q*zw - r*ywp + Vp_w;
vwp_ycomp = v + r*xw - p*zw;
vwp_zcomp = w + p*ywp - q*xw + Vp_w_swirl;
vwp = sqrt(vwp_xcomp.^2 + vwp_ycomp.^2 + vwp_zcomp.^2);
vwp_xz = sqrt(vwp_xcomp.^2 + vwp_zcomp.^2); 

vts_xcomp = u + q*zt - r*yts + Vp_t;
vts_ycomp = v + r*xt - p*zt;
vts_zcomp = w + p*yts - q*xt - Vp_t_swirl;
vts = sqrt(vts_xcomp.^2 + vts_ycomp.^2 + vts_zcomp.^2);
vts_xz = sqrt(vts_xcomp.^2 + vts_zcomp.^2);

vtp_xcomp = u + q*zt - r*ytp + Vp_t;
vtp_ycomp = v + r*xt - p*zt;
vtp_zcomp = w + p*ytp - q*xt + Vp_t_swirl;
vtp = sqrt(vtp_xcomp.^2 + vtp_ycomp.^2 + vtp_zcomp.^2);
vtp_xz = sqrt(vtp_xcomp.^2 + vtp_zcomp.^2);

vr_xcomp = u + q*zr - r*yr + Vp_r;
vr_ycomp = v + r*xr - p*zr + sign(zr + zCG).*Vp_r_swirl;
vr_zcomp = w + p*yr - q*xr;
vr = sqrt(vr_xcomp.^2 + vr_ycomp.^2 + vr_zcomp.^2);
vr_xy = sqrt(vr_xcomp.^2 + vr_ycomp.^2);

vB_xcomp = u + q*zB - r*yB + Vp_B;
vB_ycomp = v + r*xB - p*zB + sign(zB + zCG).*Vp_B_swirl;
vB_zcomp = w + p*yB - q*xB;
vB = sqrt(vB_xcomp.^2 + vB_ycomp.^2 + vB_zcomp.^2);
vB_xy = sqrt(vB_xcomp.^2 + vB_ycomp.^2);

%--------------------------------------------------------------------------
%              Wing, Tail, Rudder and Body angles of attack
%--------------------------------------------------------------------------
% Angle of attack = atan(V_zcomp / V_xcomp)
% range is -180 -> 180
a_ws = atan2(vws_zcomp,vws_xcomp);
a_wp = atan2(vwp_zcomp,vwp_xcomp);

a_ts = atan2(vts_zcomp,vts_xcomp);
a_tp = atan2(vtp_zcomp,vtp_xcomp);

% Vertical AoA for rudder and body = atan(V_ycomp / V_xcomp)
% range is -180 -> 180
a_r = atan2(vr_ycomp,vr_xcomp);
a_B  = atan2(vB_ycomp,vB_xcomp);
    
%--------------------------------------------------------------------------
%                   Wing lift and drag coefficient
%-------------------------------------------------------------------------- 
CN_ws = zeros(1,Nw);
CL_ws = zeros(1,Nw);
CD_ws = zeros(1,Nw);
CM_ws = zeros(1,Nw);

CN_wp = zeros(1,Nw);
CL_wp = zeros(1,Nw);
CD_wp = zeros(1,Nw);
CM_wp = zeros(1,Nw);

for i = 1:Nw
    [CN_ws(i),CL_ws(i),CD_ws(i),CM_ws(i)] = flappedAirfoil_LAR_v1_FlatPlate(i,a_ws(i),Cwf(i),Cw(i),RAilDef,AR_w,Cd0,Cd90,alp0);
    [CN_wp(i),CL_wp(i),CD_wp(i),CM_wp(i)] = flappedAirfoil_LAR_v1_FlatPlate(i,a_wp(i),Cwf(i),Cw(i),-RAilDef,AR_w,Cd0,Cd90,alp0);    
end

%--------------------------------------------------------------------------
%                       Wing-Tail interference
%--------------------------------------------------------------------------
% % Bound vortices of wing
% G_ws = 0.5*Cw.*CL_ws.*vws_xz;
% G_wp = 0.5*Cw.*CL_wp.*vwp_xz;
% 
% % Induced downwash on tail due to bound vortices
% ViTS_ws = zeros(Nt,Nw);
% ViTS_wp = zeros(Nt,Nw);
% ViTP_ws = zeros(Nt,Nw);
% ViTP_wp = zeros(Nt,Nw);
% 
% for i = 1:Nt
%     for j = 1:Nw
%         ViTS_ws(i,j) = G_ws(j)*GTermTS_ws(i,j);
%         ViTS_wp(i,j) = G_wp(j)*GTermTS_wp(i,j);
%         ViTP_ws(i,j) = G_ws(j)*GTermTP_ws(i,j);
%         ViTP_wp(i,j) = G_wp(j)*GTermTP_wp(i,j);    
%     end
% end
% 
% % Trailing vortices of wing
% delG_ws = zeros(1,Nw+2);
% delG_wp = zeros(1,Nw+2);
% 
% delG_ws(1) = G_ws(1) - 0;
% delG_wp(1) = 0 - G_wp(1);
% 
% for i = 2:Nw-1
%     delG_ws(i) = G_ws(i) - G_ws(i-1);
%     delG_wp(i) = G_wp(i-1) - G_wp(i);
% end
% 
% delG_ws(Nw) = G_ws(Nw);
% delG_ws(Nw+1) = 0 - G_ws(Nw-1);
% delG_ws(Nw+2) = 0 - G_ws(Nw);
% 
% delG_wp(Nw) = 0 - G_wp(Nw);
% delG_wp(Nw+1) = G_wp(Nw-1) - 0;
% delG_wp(Nw+2) = G_wp(Nw) - 0;
% 
% % Induced downwash on tail due to trailing vortices
% delViTS_ws = zeros(Nt,Nw+2);
% delViTS_wp = zeros(Nt,Nw+2);
% delViTP_ws = zeros(Nt,Nw+2);
% delViTP_wp = zeros(Nt,Nw+2);
% 
% for i = 1:Nt
%     for j = 1:Nw+2
%         delViTS_ws(i,j) = delG_ws(j)*delGTermTS_ws(i,j);
%         delViTS_wp(i,j) = delG_wp(j)*delGTermTS_wp(i,j);
%         delViTP_ws(i,j) = delG_ws(j)*delGTermTP_ws(i,j);
%         delViTP_wp(i,j) = delG_wp(j)*delGTermTP_wp(i,j);
%     end
% end
%  
% % Summing all induced downwash
% ViTS = zeros(1,Nt);
% ViTP = zeros(1,Nt);
% 
% for i = 1:Nt
%     ViTS(i) = sum(ViTS_ws(i,:)) + sum(ViTS_wp(i,:)) + sum(delViTS_ws(i,:)) + sum(delViTS_wp(i,:));
%     ViTP(i) = sum(ViTP_ws(i,:)) + sum(ViTP_wp(i,:)) + sum(delViTP_ws(i,:)) + sum(delViTP_wp(i,:));   
% end
%  
% % Induced angle of attack on tail
% vts_zcomp = vts_zcomp - ViTS;
% vtp_zcomp = vtp_zcomp - ViTP;
% 
% a_ts = atan2(vts_zcomp,vts_xcomp);
% a_tp = atan2(vtp_zcomp,vtp_xcomp);

%--------------------------------------------------------------------------
%                  Tail's lift and drag coefficient
%--------------------------------------------------------------------------
CN_ts = zeros(1,Nt);
CL_ts = zeros(1,Nt);
CD_ts = zeros(1,Nt);
CM_ts = zeros(1,Nt);

CN_tp = zeros(1,Nt);
CL_tp = zeros(1,Nt);
CD_tp = zeros(1,Nt);
CM_tp = zeros(1,Nt);

for i = 1:Nt    
    [CN_ts(i),CL_ts(i),CD_ts(i),CM_ts(i)] = flappedAirfoil_LAR_v1_FlatPlate(i,a_ts(i),Ctf(i),Ct(i),ElevDef,AR_t,Cd0,Cd90,alp0);
    [CN_tp(i),CL_tp(i),CD_tp(i),CM_tp(i)] = flappedAirfoil_LAR_v1_FlatPlate(i,a_tp(i),Ctf(i),Ct(i),ElevDef,AR_t,Cd0,Cd90,alp0);
end

%--------------------------------------------------------------------------
%                       Wing-Rudder interference
%--------------------------------------------------------------------------

% % Induced downwash on rudder due to bound vortices
% ViR_ws_xcomp = zeros(Nr,Nw);
% ViR_ws_zcomp = zeros(Nr,Nw);
% ViR_wp_xcomp = zeros(Nr,Nw);
% ViR_wp_zcomp = zeros(Nr,Nw);
% 
% for i = 1:Nr
%     for j = 1:Nw
%         ViR_ws_xcomp(i,j) = G_ws(j)*GTermR_w_xcomp(i,j);
%         ViR_ws_zcomp(i,j) = G_ws(j)*GTermR_w_zcomp(i,j);
%         
%         ViR_wp_xcomp(i,j) = G_wp(j)*GTermR_w_xcomp(i,j);
%         ViR_wp_zcomp(i,j) = G_wp(j)*GTermR_w_zcomp(i,j);    
%     end
% end
% 
% % Induced downwash on rudder due to trailing vortices
% delViR_ws_ycomp = zeros(Nr,Nw+1);
% delViR_ws_zcomp = zeros(Nr,Nw+1);
% delViR_wp_ycomp = zeros(Nr,Nw+1);
% delViR_wp_zcomp = zeros(Nr,Nw+1);
% 
% for i = 1:Nr
%     for j = 1:Nw+1
%         delViR_ws_ycomp(i,j) = delG_ws(j)*delGTermR_ws_ycomp(i,j);
%         delViR_ws_zcomp(i,j) = delG_ws(j)*delGTermR_ws_zcomp(i,j);
%         
%         delViR_wp_ycomp(i,j) = delG_wp(j)*delGTermR_wp_ycomp(i,j);        
%         delViR_wp_zcomp(i,j) = delG_wp(j)*delGTermR_wp_zcomp(i,j);
%     end
% end
% 
% % Summing all induced downwash
% ViR_xcomp = zeros(1,Nr);
% ViR_ycomp = zeros(1,Nr);
% ViR_zcomp = zeros(1,Nr);
% 
% for i = 1:Nr
%     ViR_xcomp(i) = sum(ViR_ws_xcomp(i,:)) + sum(ViR_wp_xcomp(i,:));
%     ViR_ycomp(i) = sum(delViR_ws_ycomp(i,:)) + sum(delViR_wp_ycomp(i,:));
%     ViR_zcomp(i) = sum(ViR_ws_zcomp(i,:)) + sum(ViR_wp_zcomp(i,:)) + sum(delViR_ws_zcomp(i,:)) + sum(delViR_wp_zcomp(i,:)); 
% end
% 
% % Induced angle of attack on rudder
% vr_xcomp = vr_xcomp - ViR_xcomp;
% vr_ycomp = vr_ycomp - ViR_ycomp;
% vr_zcomp = vr_zcomp - ViR_zcomp;
% 
% a_r = atan2(vr_ycomp,vr_xcomp);

%--------------------------------------------------------------------------
%                  Rudder lift and drag coefficient
%--------------------------------------------------------------------------
CN_r = zeros(1,Nr);
CL_r = zeros(1,Nr);
CD_r = zeros(1,Nr);
CM_r = zeros(1,Nr);

for i = 1:Nr    % check sign of deflection into flapped code
    [CN_r(i),CL_r(i),CD_r(i),CM_r(i)] = flappedAirfoil_LAR_v1_FlatPlate(i,a_r(i),Crf(i),Cr(i),-RudDef,AR_r,Cd0,Cd90,alp0);
end

%--------------------------------------------------------------------------
%                  Fuselage lift and drag coefficient
%--------------------------------------------------------------------------
CN_B = zeros(1,NB);
CL_B = zeros(1,NB);
CD_B = zeros(1,NB);
CM_B = zeros(1,NB);

for i = 1:NB
     [CN_B(i),CL_B(i),CD_B(i),CM_B(i)] = flappedAirfoil_LAR_v1_FlatPlate(i,a_B(i),CBf(i),CB(i),0,AR_B,Cd0,Cd90,alp0);
end

%--------------------------------------------------------------------------
%                   Aerodynamic Forces & Moments
%--------------------------------------------------------------------------
% x-direction forces
Fx_ws = 0.5*rho*bw.*Cw.*vws_xz.^2.*(CL_ws.*sin(a_ws) - CD_ws.*cos(a_ws));
Fx_wp = 0.5*rho*bw.*Cw.*vwp_xz.^2.*(CL_wp.*sin(a_wp) - CD_wp.*cos(a_wp));
Fx_ts = 0.5*rho*bt.*Ct.*vts_xz.^2.*(CL_ts.*sin(a_ts) - CD_ts.*cos(a_ts));
Fx_tp = 0.5*rho*bt.*Ct.*vtp_xz.^2.*(CL_tp.*sin(a_tp) - CD_tp.*cos(a_tp));

Fx_r = 0.5*rho*br.*Cr.*vr_xy.^2.*(CL_r.*sin(a_r) - CD_r.*cos(a_r));
Fx_B = 0.5*rho*bB.*CB.*vB_xy.^2.*(CL_B.*sin(a_B) - CD_B.*cos(a_B));

% y-direction forces
Fy_ws = [0,0,0,0,0,0,0]; % INCLUDE FRICTION DRAG
Fy_wp = [0,0,0,0,0,0,0]; % INCLUDE FRICTION DRAG
Fy_ts = [0,0,0]; % INCLUDE FRICTION DRAG
Fy_tp = [0,0,0]; % INCLUDE FRICTION DRAG

Fy_r = 0.5*rho*br.*Cr.*vr_xy.^2.*(-CL_r.*cos(a_r) - CD_r.*sin(a_r));
Fy_B = 0.5*rho*bB.*CB.*vB_xy.^2.*(-CL_B.*cos(a_B) - CD_B.*sin(a_B));

% z-direction forces
Fz_ws = 0.5*rho*bw.*Cw.*vws_xz.^2.*(-CL_ws.*cos(a_ws) - CD_ws.*sin(a_ws));
Fz_wp = 0.5*rho*bw.*Cw.*vwp_xz.^2.*(-CL_wp.*cos(a_wp) - CD_wp.*sin(a_wp));
Fz_ts = 0.5*rho*bt.*Ct.*vts_xz.^2.*(-CL_ts.*cos(a_ts) - CD_ts.*sin(a_ts));
Fz_tp = 0.5*rho*bt.*Ct.*vtp_xz.^2.*(-CL_tp.*cos(a_tp) - CD_tp.*sin(a_tp));
Fz_r = [0,0,0,0]; % INCLUDE FRICTION DRAG
Fz_B = [0,0,0,0]; % INCLUDE FRICTION DRAG

% x-direction moments
Mx_ws = yws.*Fz_ws - zw.*Fy_ws;
Mx_wp = ywp.*Fz_wp - zw.*Fy_wp;
Mx_ts = yts.*Fz_ts - zt.*Fy_ts;
Mx_tp = ytp.*Fz_tp - zt.*Fy_tp;
Mx_r = yr.*Fz_r - zr.*Fy_r;
Mx_B = yB.*Fz_B - zB.*Fy_B;
Mx_Thr = yp*ThrForce(3) - zp*ThrForce(2);

% y-direction moments
My_ws = zw.*Fx_ws - xw.*Fz_ws + 0.5*rho*bw.*Cw.*vws_xz.^2.*Cw.*CM_ws;
My_wp = zw.*Fx_wp - xw.*Fz_wp + 0.5*rho*bw.*Cw.*vwp_xz.^2.*Cw.*CM_wp;
My_ts = zt.*Fx_ts - xt.*Fz_ts + 0.5*rho*bt.*Ct.*vts_xz.^2.*Ct.*CM_ts;
My_tp = zt.*Fx_tp - xt.*Fz_tp + 0.5*rho*bt.*Ct.*vtp_xz.^2.*Ct.*CM_tp;
My_r = zr.*Fx_r - xr.*Fz_r;
My_B = zB.*Fx_B - xB.*Fz_B;
My_Thr = zp*ThrForce(1) - xp*ThrForce(3);

% z-direction moments
Mz_ws = xw.*Fy_ws - yws.*Fx_ws;
Mz_wp = xw.*Fy_wp - ywp.*Fx_wp;
Mz_ts = xt.*Fy_ts - yts.*Fx_ts;
Mz_tp = xt.*Fy_tp - ytp.*Fx_tp;
Mz_r = xr.*Fy_r - yr.*Fx_r - 0.5*rho*br.*Cr.*vr_xy.^2.*Cr.*CM_r;  % WHY NEGATIVE SIGN???
Mz_B = xB.*Fy_B - yB.*Fx_B - 0.5*rho*bB.*CB.*vB_xy.^2.*CB.*CM_B;  % WHY NEGATIVE SIGN???
Mz_Thr = xp*ThrForce(2) - yp*ThrForce(1);

% Total forces and moments in x, y and z directions
Fx = sum(Fx_ws) + sum(Fx_wp) + sum(Fx_ts) + sum(Fx_tp) + sum(Fx_r) + sum(Fx_B) + ThrForce(1);
Fy = sum(Fy_ws) + sum(Fy_wp) + sum(Fy_ts) + sum(Fy_tp) + sum(Fy_r) + sum(Fy_B) + ThrForce(2);
Fz = sum(Fz_ws) + sum(Fz_wp) + sum(Fz_ts) + sum(Fz_tp) + sum(Fz_r) + sum(Fz_B) + ThrForce(3);

Mx = sum(Mx_ws) + sum(Mx_wp) + sum(Mx_ts) + sum(Mx_tp) + sum(Mx_r) + sum(Mx_B) + Mx_Thr + ThrMom(1) + GyroMom(1);
My = sum(My_ws) + sum(My_wp) + sum(My_ts) + sum(My_tp) + sum(My_r) + sum(My_B) + My_Thr + ThrMom(2) + GyroMom(2);
Mz = sum(Mz_ws) + sum(Mz_wp) + sum(Mz_ts) + sum(Mz_tp) + sum(Mz_r) + sum(Mz_B) + Mz_Thr + ThrMom(3) + GyroMom(3);

a= 0;