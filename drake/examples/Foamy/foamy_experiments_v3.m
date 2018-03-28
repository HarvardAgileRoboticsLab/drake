clear all; close all; clc;

oldpath = path;
path = path(oldpath, 'Controllers');

% create a plane model and run a simple trajectory (trim conditions)
plane = FoamyPlant();
[xtraj, utraj] = runDircol(plane, 0);

tspan = xtraj(1).tspan;
t = linspace(tspan(1), tspan(2), 100);

% try running TVLQR from a different initial condition towards the desired trajectory
Q = 50 * eye(12);
R = eye(4);
Qf = Q;
t0 = xtraj.tspan(1);
tf = xtraj.tspan(2);

% initial conditions with disturbance
x00 = xtraj.eval(0);
x00(1:3) = x00(1:3) + [0.05; 0.05; -0.05];

% build all controllers
control_discrete = FoamyDiscreteTVLQRController(plane, xtraj, utraj, Q, R, Qf, t0, tf);
control_tvlqr = FoamyTVLQRController(plane, xtraj, utraj, Q, R, Qf, t0, tf);
control_phi = FoamyPhiTVLQRController(plane, xtraj, utraj, Q, R, Qf, t0, tf);

% create feedback models
cplane_discrete = ControlledFoamyPlant(control_discrete);
cplane_tvlqr = ControlledFoamyPlant(control_tvlqr);
cplane_phi = ControlledFoamyPlant(control_phi);

% simulate all systems
tf_phi = 2.2;
xtraj_discrete = simulate(cplane_discrete, [t0 tf], x00);
xtraj_tvlqr = simulate(cplane_tvlqr, [t0 tf], x00);
xtraj_phi = simulate(cplane_phi, [t0 tf_phi], x00);

% plot trajectories
plotting(t0, tf, tf_phi, xtraj, xtraj_discrete, xtraj_tvlqr, xtraj_phi);


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
Qf2 = 10 * Q; %Qf2(1:3, 1:3) = 10 * Qf2(1:3, 1:3); % weigh the position highly

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


%% plot K and S matrices across controllers

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
  
% plot phi-tvlqr K values
% subplot(2, 3, 3),
% plot([control_phi_1.K_array.t, control_phi_2.K_array.t],...
%      [control_phi_1.K_array.K, control_phi_2.K_array.K]), hold on
% plot([control_phi_1.tf, control_phi_1.tf],...
%      [min(control_phi_2.K_array.K(:)), 1.5 * max(control_phi_2.K_array.K(:))], 'r--'), hold off
% title('TVLQR: K');

 % plot phi-tvlqr S values
% subplot(2, 3, 6)
% plot([control_phi_1.S_array.t, control_phi_2.S_array.t],...
%      [control_phi_1.S_array.S, control_phi_2.S_array.S]), hold on
% plot([control_phi_1.tf, control_phi_1.tf],...
%      [min(control_phi_2.S_array.S(:)), 1.5 * max(control_phi_2.S_array.S(:))], 'r--'), hold off
% title('TVLQR: S');

% figure of matrix eigenvalues
figure
subplot(1, 3, 1), hold on
for i = 1:numel(control_discrete_2.S_array.t)
    plot(control_discrete_2.S_array.t(i) * ones(12, 1), eig(reshape(control_discrete_2.S_array.S(:, i), [12 12])), 'r.');
end
for i = 1:numel(control_tvlqr_2.S_array.t)
    plot(control_tvlqr_2.S_array.t(i) * ones(12, 1), eig(reshape(control_tvlqr_2.S_array.S(:, i), [12 12])), 'b.');
end

%% -- phase-based two-part controller --

% control_phi_1 = FoamyPhiTVLQRController(...
%     hplane.modes{1}, xtraj.traj{1}.trajs{2}, utraj.traj{1}, Q, R, Qf1, t0, xtraj.te);
% 
% ch_phi = ControlledHybridFoamyPlant(control_phi_1, control_tvlqr_2);


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

%%
ABtest = ABTestSystem(hplane.modes{1}, utraj.traj{1});

            
% do quaternion transformation
x1 = xtraj.traj{1}.trajs{2}.eval(0); x1 = [x1(1:3); x1(5:7); x1(8:13)];
xx = simulate(ABtest, [t0 xtraj.te], x1);



