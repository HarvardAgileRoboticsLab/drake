clear; clc; close all; 
global kl_traj u_traj 
kl_traj = []; 
u_traj = []; 

%% Load Trajectory

load TrajOpt-FixedBody_SimpleSprings_fullRobot

tt = xtraj.FR_scaled.getBreaks(); 
hh = mean(diff(tt)); 
uu_des = utraj.FR_scaled.eval(tt + hh/2); 

%% Build SL
sl_urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'dev', 'SimpleHAMR', 'urdf', ...
    'FR_scaled.urdf');

% options
options.terrain = [];
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = false;
options.pf = [0, -7.58, -11.35]';
options.foot = 'FRL4';
options.dt = hh; 

SL = SLTSRBM(sl_urdf, options);
v = SL.constructVisualizer();

nq = SL.getNumPositions();
nv = SL.getNumVelocities();
nx = nq + nv;
nu = SL.getNumInputs();
nl = SL.nl;

%% Simulate

utraj_sim = PPTrajectory(zoh(tt, uu_des(1:nu, :))); 
utraj_sim = utraj_sim.setOutputFrame(SL.getInputFrame()); 
SL_OL = cascade(utraj_sim, SL); 
xtraj_sim = simulate(SL_OL, [tt(1), tt(end)], xtraj.FR_scaled.eval(tt(1))); 

%%
ll_des = uu_des(nu+(1:nl), :);
xx_des = xtraj.FR_scaled.eval(tt); 
xx_sim = xtraj_sim.eval(tt); 

act_dof = SL.getActuatedJoints();

figure(1); clf; hold on; 
for i = 1:nu
    subplot(2,1,i); hold on; 
    plot(tt, u_traj(i,:)); 
    plot(tt + hh/2, uu_des(i,:))    
    legend('Input', 'Desired Input')
end

figure(2); clf; hold on;
for i = 1:numel(act_dof)
    subplot(2,1,i); hold on;
    plot(tt + hh, xx_sim(act_dof(i), :)*1e3);
    plot(tt, xx_des(act_dof(i), :)*1e3);
    legend('Act Deflection', 'Desired Act. Deflection');
end

figure(3); clf; hold on;
for i = 1:size(ll_des, 1)
    subplot(2,3,i); hold on;
    plot(tt + hh, kl_traj(i, :));
    plot(tt + hh/2, ll_des(i,:));
    legend('Constraint Force', 'Desired Constraint Force');
end
