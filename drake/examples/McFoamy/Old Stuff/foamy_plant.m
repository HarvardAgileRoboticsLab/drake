function [qdot] = foamy_plant(t,x,u)
%x has 14 dimensions
Quat = x(4:7);
Vel = x(8:10);
B_rate = x(11:13);

ThrRPM = u(1);
Aileron = u(2);
Elevator = u(3);
Rudder = u(4);


%Initial Conditions
%global Init_Quat Init_Pos
Init_Pos=[0 0 0];
Init_Vel=[5 0 0];
Init_Eul=[0 10 0]*pi/180;
Init_Quat=EulToQuat(Init_Eul);
Init_Omega=[0 0 0];

%No wind for now
V_wind = zeros(3,1);

%Vel = Init_Vel;
%B_rate = Init_Quat;
%Quat = Init_Quat;

%Maximum Actuators
AilDef_max = 40; %degrees
ElevDef_max = 50; %degrees
RudDef_max = 50; %degrees
ThrRPM_max = 6710; %RPM

Controls = [Aileron;Elevator;Rudder;ThrRPM];

%Check if command has saturated in either direction
[AilDefOut,ElevDefOut,RudDefOut,ThrRPMOut] = Command_saturation(AilDef_max,ElevDef_max,RudDef_max,ThrRPM_max,Controls);

%Control Definitions
CtrlDef = (pi/180).*[AilDefOut;-1*AilDefOut;ElevDefOut;RudDefOut];

%Aircraft Parameters
[PropGeom,Geometry,CG] = Aircraft_parameters();

%Calculate the Thruster Forces and Moments
[ThrForce,ThrMom,RPM,J,Vi0_avg] = Thruster_forces_and_moments(ThrRPMOut, Vel, B_rate, PropGeom, CG);

%Calculate the Gyroscopic Moment
[GyroMom] = Gyro_moment(RPM*(2*pi/60),B_rate,J);

%Calculate prop velocities w/slipstream model
[VProp_Axial, VProp_Swirl] = Slipstream_model(RPM, Vi0_avg, PropGeom, Geometry);

%Aerodynamics calculations for forces and moments
[Fx,Fy,Fz,Mx,My,Mz] = Aerodynamics(Vel,CtrlDef,ThrForce,ThrMom,GyroMom,VProp_Axial,VProp_Swirl,V_wind,Geometry,CG,B_rate);

%Boddy accelerations calculation
[qdd] = Body_accelerations(Vel,Fx,Fy,Fz,Mx,My,Mz,Quat,B_rate);
qdot = [Vel;B_rate;qdd];
end

