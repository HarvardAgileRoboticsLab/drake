%% HYBRID


oldpath = path();
addpath('Controllers', 'Controllers/PendulumControllers');

clear all; close all; clc;

% make hybrid pendulum plant
mass = 0.1;
plant = HybridFoamyPendulumPlant(mass);

% % load and evaluate trajectory
load('/home/user/drake/drake/examples/Foamy/Testing/Traj_Test/test3.mat');
t0 = xtraj.tspan(1);
tf = xtraj.tspan(2);
t = linspace(t0, tf, 100);
xtraj_eval = xtraj.eval(t);

% plot trajectory
figure
plot3(xtraj_eval(2,:), xtraj_eval(3,:), xtraj_eval(4,:), 'r'); hold on;
title('Plane trajectory in space');
xlabel('x'), ylabel('y'), zlabel('z');
axis equal

% set TVLQR matrices
Q = diag([10 * ones(1, 12), 0 * ones(1, 12)]);
R = eye(4);
Qf2 = Q;

% take state array before and after collision
xx_post = xtraj.traj{2}.trajs{2};
xx_pre = xtraj.traj{1}.trajs{2};

% calculate post-collision controller
control_tvlqr_2 = PendulumTVLQRController(...
     plant.modes{2}, xx_post, utraj.traj{2}, Q, R, Qf2, xtraj.te, tf);
 
%  ================
% find impact mapping for S
dt = 0.05;
te = xtraj.te;
xtrans = xtraj.traj{1}.trajs{2}.eval(te);
utrans = utraj.traj{1}.eval(te);
[xdot,dxdot] = plant.impact_dynamics(te, xtrans, utrans, dt);
A_impact_quat = dxdot(:, 2:27);

hat = @(x) [  0   -x(3)  x(2)
         x(3)   0   -x(1)
        -x(2)  x(1)  0];

% first Gk
qk = xtrans(4:7);
sk = qk(1);
vk = qk(2:4);

Gk = [-vk'; sk*eye(3) + hat(vk)];

% construct second Gk (for pendulum)
qp = xtrans(11:14);
sp = qp(1);
vp = qp(2:4);

Gp = [-vp'; sp*eye(3) + hat(vp)];

% discrete matrices?
convert = blkdiag(eye(3), Gk, eye(3), Gp, eye(12));
A_impact = convert' * A_impact_quat * convert;
% ================


Qf1 = A_impact' * reshape(control_tvlqr_2.S_array.S(:, 1), [24 24]) * A_impact;

control_tvlqr_1 = PendulumTVLQRController(...
    plant.modes{1}, xx_pre, utraj.traj{1}, Q, R, Qf1, t0, xtraj.te);

pendulum_system = ControlledHybridFoamyPendulumPlant(control_tvlqr_1, control_tvlqr_2, mass);

%% simulate

% initial conditions with disturbance
% state: [m, x, y, z, q1, q2, q3, q4, vx, vy, vz, wx, wy, wz];
x00 = xtraj.eval(0);
x00(2:4) = x00(2:4) + [-0.2; 0.2; -0.2];
x00(9:11) = x00(9:11) + [-0.2; 0.2; -0.2];

% simulate
xtraj_tvlqr = simulate(pendulum_system, [t0 tf], x00);
xtraj_tvlqr_eval = xtraj_tvlqr.eval(t);

% check guard values
% tt_arr = 0.99:0.00001:1.0001;
% for i = 1:numel(tt_arr)
%     tt = tt_arr(i);
%     x = xtraj.eval(tt); x = x(2:end);
%     gz = plant.packageGuardz(tt, x, utraj.eval(tt));
%     gx = plant.packageGuardx(tt, x, utraj.eval(tt));
%     dist = norm(x(1:3));
%     fprintf('t = %.5f, gx = %.3f, gz = %.3f, dist = %.3f\n', tt, gx, gz, dist);
% end

% calculate switch points
xtraj_switch = xtraj.eval(xtraj.te);
xtraj_tvlqr_switch = xtraj_tvlqr.eval(xtraj_tvlqr.te);

% plot 3D positions for all controller trajectories
figure
hold on
plot3(xtraj_eval(2,:), xtraj_eval(3,:), xtraj_eval(4,:), 'r')
plot3(xtraj_tvlqr_eval(2,:), xtraj_tvlqr_eval(3,:), xtraj_tvlqr_eval(4,:), 'k')
scatter3(xtraj_switch(2,:), xtraj_switch(3,:), xtraj_switch(4,:), 'ro')
scatter3(xtraj_tvlqr_switch(2,:), xtraj_tvlqr_switch(3,:), xtraj_tvlqr_switch(4,:), 'ko')
hold off
title('Plane trajectory in space');
xlabel('x'), ylabel('y'), zlabel('z');
legend('nominal', 'TVLQR', 'Location', 'northeast');
axis equal

% plot all states for all controller trajectories
states = {'mode', 'plane x', 'plane y', 'plane z', 'plane q1', 'plane q2', 'plane q3', 'plane q4',...
                  'hook x', 'hook y', 'hook z', 'hook q1', 'hook q2', 'hook q3', 'hook q4',...
                  'plane vx', 'plane vy', 'plane vz', 'plane w1', 'plane w2', 'plane w3',...
                  'hook vx', 'hook vy', 'hook vz', 'hook w1', 'hook w2', 'hook w3'};
figure('Position', [100 100 1000 800]);

subplot(5, 7, 4)
plot(t, xtraj_eval(1,:), 'r'), hold on
plot(t, xtraj_tvlqr_eval(1,:), 'k'), hold off
title(states{1}, 'Fontsize', 14);
set(gca, 'Xlim', [t0 tf]);
ytickformat('%.2f');
for i = 2:26
    subplot(5, 7, 7 + i)
    plot_num = i;
    plot(t, xtraj_eval(plot_num,:), 'r'), hold on
    plot(t, xtraj_tvlqr_eval(plot_num,:), 'k'), hold off
    title(states{plot_num}, 'Fontsize', 14);
    set(gca, 'Xlim', [t0 tf]);
    ytickformat('%.2f');
end
