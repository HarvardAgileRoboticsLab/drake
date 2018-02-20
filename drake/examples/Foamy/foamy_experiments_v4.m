%% HYBRID

clear all; close all; clc;

% create a plane model and run a simple trajectory (trim conditions)
hplane = HybridFoamyPlant();
[xtraj, utraj] = runDircol(hplane, 0);

tspan = xtraj(1).tspan;
t = linspace(tspan(1), tspan(2), 100);

% try running TVLQR from a different initial condition towards the desired trajectory
Q = 10 * eye(12);
R = eye(4);
Qf2 = Q; Qf2(1:3, 1:3) = 10 * Qf2(1:3, 1:3); % weigh the position highly

t0 = xtraj.tspan(1);
tf = xtraj.tspan(2);


%% -- discrete tvlqr two-part controller --

control_discrete_2 = FoamyDiscreteTVLQRController(...
     hplane.modes{2}, xtraj.traj{2}.trajs{2}, utraj.traj{2}, Q, R, Qf2, xtraj.te, tf);
 
Qf1 = reshape(control_discrete_2.S_array.S(:, 1), [12 12]);

control_discrete_1 = FoamyDiscreteTVLQRController(...
     hplane.modes{1}, xtraj.traj{1}.trajs{2}, utraj.traj{1}, Q, R, Qf1, t0, xtraj.te);
 
ch_discrete = ControlledHybridFoamyPlant(control_discrete_1, control_discrete_2);


%% -- tvlqr two-part controller --
control_tvlqr_2 = FoamyTVLQRController(...
     hplane.modes{2}, xtraj.traj{2}.trajs{2}, utraj.traj{2}, Q, R, Qf2, xtraj.te, tf);

Qf1 = reshape(control_tvlqr_2.S_array.S(:, 1), [12 12]);

control_tvlqr_1 = FoamyTVLQRController(...
    hplane.modes{1}, xtraj.traj{1}.trajs{2}, utraj.traj{1}, Q, R, Qf1, t0, xtraj.te);

ch_tvlqr = ControlledHybridFoamyPlant(control_tvlqr_1, control_tvlqr_2);


%% plotting time!

figure

% plot discrete K values
subplot(2, 2, 1),
plot([control_discrete_1.K_array.t, control_discrete_2.K_array.t],...
     [control_discrete_1.K_array.K, control_discrete_2.K_array.K]), hold on
plot([control_discrete_1.tf, control_discrete_1.tf],...
     [min(control_discrete_2.K_array.K(:)), 1.5 * max(control_discrete_2.K_array.K(:))], 'r--'), hold off
title('Discrete TVLQR: K');
 
 % plot discrete S values
subplot(2, 2, 2)
plot([control_discrete_1.S_array.t, control_discrete_2.S_array.t],...
     [control_discrete_1.S_array.S, control_discrete_2.S_array.S]), hold on
plot([control_discrete_1.tf, control_discrete_1.tf],...
     [min(control_discrete_2.S_array.S(:)), 1.5 * max(control_discrete_2.S_array.S(:))], 'r--'), hold off
title('Discrete TVLQR: S');
 
% plot my tvlqr K values
subplot(2, 2, 3),
plot([control_tvlqr_1.K_array.t, control_tvlqr_2.K_array.t],...
     [control_tvlqr_1.K_array.K, control_tvlqr_2.K_array.K]), hold on
plot([control_tvlqr_1.tf, control_tvlqr_1.tf],...
     [min(control_tvlqr_2.K_array.K(:)), 1.5 * max(control_tvlqr_2.K_array.K(:))], 'r--'), hold off
title('TVLQR: K');
 
 % plot my tvlqr S values
subplot(2, 2, 4)
plot([control_tvlqr_1.S_array.t, control_tvlqr_2.S_array.t],...
     [control_tvlqr_1.S_array.S, control_tvlqr_2.S_array.S]), hold on
plot([control_tvlqr_1.tf, control_tvlqr_1.tf],...
     [min(control_tvlqr_2.S_array.S(:)), 1.5 * max(control_tvlqr_2.S_array.S(:))], 'r--'), hold off
title('TVLQR: S');


%% simulate

% initial conditions with disturbance
x00 = xtraj.eval(0);
% disturbance: [m, x, y, z, q1, q2, q3, q4, vx, vy, vz, wx, wy, wz];
x00 = x00 + [0; -1; 1; -1; 0; 0; 0; 0; 0.1; 0.1; 0; 0; 0; 0];

% run all simulations
%xtraj_phi = simulate(ch_phi, [t0 tf], x00);
xtraj_tvlqr = simulate(ch_tvlqr, [t0 tf], x00);
xtraj_discrete = simulate(ch_discrete, [t0 tf], x00);

xtraj_phi = xtraj_tvlqr;
%xtraj_discrete = xtraj_tvlqr;

% plot trajectories
plottingHybrid(t0, tf, tf, xtraj, xtraj_discrete, xtraj_tvlqr, xtraj_phi);



%% PENDULUM VERSION 

clear all; close all; clc;

% load and evaluate trajectory
load('/home/user/drake/drake/examples/Foamy/Testing/Traj_Test/test3.mat');
t = linspace(xtraj.tspan(1), xtraj.tspan(2), 100);
xtraj_eval = xtraj.eval(t);

% plot trajectory
figure
plot3(xtraj_eval(2,:), xtraj_eval(3,:), xtraj_eval(4,:), 'r'); hold on;
title('Plane trajectory in space');
xlabel('x'), ylabel('y'), zlabel('z');
axis equal

% make hybrid pendulum plant
mass = 0.1;
plant = HybridFoamyPendulumPlant(mass);

t0 = xtraj.tspan(1);
tf = xtraj.tspan(2);
t = linspace(t0, tf, 100);

% set TVLQR matrices
Q = 10 * eye(12);
R = eye(4);
Qf2 = Q; Qf2(1:3, 1:3) = 10 * Qf2(1:3, 1:3); % weigh the position highly

% take non-pendulum part of state
xx_post = xtraj.traj{2}.trajs{2}; %xx_post = xx_post(1:13);
xx_pre = xtraj.traj{1}.trajs{2}; %xx_pre = xx_pre(1:13);

% create pre-collision and post-collision controllers
control_tvlqr_2 = FoamyTVLQRController(...
     plant.modes{2}, xx_post, utraj.traj{2}, Q, R, Qf2, xtraj.te, tf);

Qf1 = reshape(control_tvlqr_2.S_array.S(:, 1), [12 12]);

control_tvlqr_1 = FoamyTVLQRController(...
    plant.modes{1}, xx_pre, utraj.traj{1}, Q, R, Qf1, t0, xtraj.te);

ch_tvlqr = ControlledHybridFoamyPlant(control_tvlqr_1, control_tvlqr_2);


% create controlled version of hybrid plant
control_hybrid_plant = ControlledHybridFoamyPendulumPlant(mass);


