clear all; close all; clc;

oldpath = path;
%path = path(oldpath, 'Controllers');
path = path(oldpath, 'Controllers/NonPendulumControllers');

% create a plane model and run a simple trajectory (trim conditions)
plane = FoamyPlant();
[xtraj, utraj] = runDircol(plane, 0);

tspan = xtraj(1).tspan;
t = linspace(tspan(1), tspan(2), 100);

% try running TVLQR from a different initial condition towards the desired trajectory
Q = 10 * eye(12);
R = eye(4);
Qn = Q;
t0 = xtraj.tspan(1);
tf = xtraj.tspan(2);

% initial conditions with disturbance
x00 = xtraj.eval(0);
x00(1:3) = x00(1:3) + [0.05; 0.05; -0.05];

% simulate discrete implementation
control = FoamyDiscreteTVLQRController(plane, xtraj, utraj, Q, R, Qn, t0, tf);
discrete_sys = feedback(plane, control);
xtraj_discrete = simulate(discrete_sys, [t0 tf], x00);

% simulate my implementation
control = FoamyTVLQRController(plane, xtraj, utraj, Q, R, Qn, t0, tf);
tvlqr_sys = feedback(plane, control);
xtraj_tvlqr = simulate(tvlqr_sys, [t0 tf], x00);

% now try phase-based implementation
%control_phi = FoamyPhiTVLQRController(plane, xtraj, utraj, Q, R, Qn, t0, tf);
%phi_sys = feedback(plane, control_phi);
%tf_phi = 2.2;
%xtraj_phi = simulate(phi_sys, [t0 tf_phi], x00);

% plot trajectories
plotting(t0, tf, tf_phi, xtraj, xtraj_discrete, xtraj_tvlqr, xtraj_tvlqr);


%% HYBRID
clear all; close all; clc;

% create a plane model and run a simple trajectory (trim conditions)
plane = HybridFoamyPlant();
[xtraj, utraj] = runDircol(plane, 0);

tspan = xtraj(1).tspan;
t = linspace(tspan(1), tspan(2), 100);

% try running TVLQR from a different initial condition towards the desired trajectory
Q = eye(12);
R = eye(4);
Qn = Q;
t0 = xtraj.tspan(1);
tf = xtraj.tspan(2);

% initial conditions with disturbance
x00 = xtraj.eval(0);
x00(2:4) = x00(2:4) + [0.05; 0.001; 0.001];

cplane = HybridControllerFeedback(control, control);



% simulate
control = HybridFoamyTVLQRController(plane, xtraj, utraj, Q, R, Qn, t0, tf);
tvlqr_sys = feedback(plane, control);
%tvlqr_sys = HybridFoamyTVLQRController(plane, xtraj, utraj, Q, R, Qn, t0, tf);
xtraj_tvlqr = simulate(tvlqr_sys, [t0 tf], x00);
xtraj_phi = xtraj_tvlqr;

% try a hack
cplane = HybridControllerFeedback(control, utraj);
xx = simulate(cplane, [t0 tf], xtraj.eval(0));

%
cplane = HybridControllerFeedback(control, control);

% evaluate trajectories
tphi = linspace(t0, 2.2, 100);
xtraj_tvlqr_eval = xtraj_tvlqr.eval(t);
xtraj_phi_eval = xtraj_phi.eval(t);
xtraj_eval = xtraj.eval(t);

% plot
states = {'mode', 'x', 'y', 'z', 'q1', 'q2', 'q3', 'q4', 'vx', 'vy', 'vz', 'w1', 'w2', 'w3'};


figure('Position', [100 100 1000 800]);
for i = 1:13
    subplot(4, 4, i)
    plot(t, xtraj_eval(i,:), 'r'), hold on
    plot(t, xtraj_tvlqr_eval(i,:), 'k'),
    plot(t, xtraj_phi_eval(i,:), 'b'), hold off
    title(states{i}, 'Fontsize', 14);
    set(gca, 'Xlim', [t0 2.2]);
    ytickformat('%.2f');
end
subplot(4, 4, [15 16])
plot(t, xtraj_eval(14,:), 'r'), hold on
plot(t, xtraj_tvlqr_eval(14,:), 'k'),
plot(t, xtraj_phi_eval(14,:), 'b'), hold off
title(states{13}, 'Fontsize', 14);
set(gca, 'Xlim', [t0 2.2]);
legend('nominal', 'TVLQR', 'phi-TVLQR', 'Location', 'bestoutside');
ytickformat('%.2f');

figure
plot3(xtraj_eval(2,:), xtraj_eval(3,:), xtraj_eval(4,:), 'r'); hold on;
plot3(xtraj_tvlqr_eval(2,:), xtraj_tvlqr_eval(3,:), xtraj_tvlqr_eval(4,:), 'k');
plot3(xtraj_phi_eval(2,:), xtraj_phi_eval(3,:), xtraj_phi_eval(4,:), 'b'); hold off;
set(gca, 'xlim', [-6 6], 'ylim', [-6 6], 'zlim', [-6 6]);
title('Plane trajectory in space');
xlabel('x'), ylabel('y'), zlabel('z')
