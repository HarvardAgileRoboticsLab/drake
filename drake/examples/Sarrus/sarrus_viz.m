clear; clc; close all;

%% Load Rigid Body

urdf = 'Sarrus.urdf'; 

% options
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.z_inactive_guess_tol = 1;
options.use_bullet = false;  
ISFLOAT = false;

% Build robot + visualizer
sarrus = RigidBodyManipulator(urdf); 
sarrus = sarrus.setGravity([0; 0; 9.81e-3]);
sarrus = compile(sarrus); 
tsarrus = TimeSteppingRigidBodyManipulator(sarrus, 1, options); 
% sarrus.manip = sarrus.manip.setGravity(sarrus.grav);

% figure(100); 
% plot(0:1:2000, -0.015 + 0.015*sin(1e-3*pi*t)


x0 = zeros(sarrus.getNumStates(), 1);

% [tf, err_str] = valuecheck(positionConstraints(sarrus,x0),zeros(18,1),1e-8);

b0 = pi/4; 
t0 = 2*b0;
o0 = -b0; 
bot = [-b0; b0; b0; -b0]; 
top = [t0; -t0; -t0; t0]; 
x0(1:sarrus.getNumDOF) = [bot; top; o0]; 

% % 
xtraj = tsarrus.simulate([0 2200], x0); 
xtraj_scaled = DTTrajectory(xtraj.getBreaks()*1e-3, xtraj.eval(xtraj.getBreaks())); 
xtraj_scaled = xtraj_scaled.setOutputFrame(xtraj.getOutputFrame()); 
v = tsarrus.constructVisualizer();
% v.inspector(x0); 

poptions.slider = true;
v.playback(xtraj_scaled, poptions);
%%
tt = xtraj.getBreaks();
yy = xtraj.eval(tt); 
q = yy(1:sarrus.getNumDOF,:); 
v = yy(sarrus.getNumDOF+1:sarrus.getNumStates, :); 
output_pos = zeros(3, numel(tt)); 
output_id = sarrus.findLinkId('Output'); 

for i = 1:numel(tt) 
    kinsol = doKinematics(sarrus, q(:,i), v(:,i)); 
    output_pos(:,i) = forwardKin(sarrus, kinsol, output_id , [2.5; 0; 0]);
end
z = output_pos(3,:);

figure(1); clf; hold on; 
subplot(2,1,1); hold on; 
title('Abs. Value of Change in Joint Angle')
plot(tt, rad2deg(abs(q(1,:)) - abs(q(1,end))))
plot(tt, rad2deg(abs(q(5,:)) - abs(q(5,end))))
plot(tt, rad2deg(abs(q(9,:)) - abs(q(9,end))), '--')
subplot(2,1,2); hold on; 
title('Height')
plot(tt, z)
% z(end) - z(1)

% 
% f = -0.015 - 0.015*sin(1e-3*pi*tt);
% f(tt <= 1000) = 0; 
% 
% figure(100); 
% plot(tt, f); 
% 
% dz = abs(z(1) - z);
% k = (0.015*tt/1000)./dz; 

% figure(1); clf; hold on; 
% subplot(1,2,1); 
% yyaxis left; plot(tt, dz)
% yyaxis right; plot(tt, f)

% 
% subplot(1,2,2); 
% plot(tt(1001:end), k(1001:end))


% plot(