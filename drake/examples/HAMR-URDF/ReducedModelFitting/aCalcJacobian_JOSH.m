clear; clc; close all;
addpath('../')
% save_dir = '~/Dropbox/CurrentWork/FrictionTrajOpt/MatFiles/SysIDFiles/';

global u_traj kl_traj
u_traj = []; 
kl_traj = []; 

% Sample Rate and Time
dt = 0.25;
Fs = 1/dt;
T = 2000;
t = 0:dt:T;
N = numel(t);

% White Noise Freq Characteristics
f1 = 40e-3;
band = [0, 2*f1*dt];

%% Build robot

% options
name = 'FL_scaled';
urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'urdf', [name, '.urdf']);
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = false;
options.z_inactive_guess_tol = 0.1;
options.dt = dt;

SL = SLTSRBM(urdf, options);
nq = SL.getNumPositions();
nv = SL.getNumVelocities();
nu = SL.getNumInputs();

% Relax joint limits
SL = SL.setJointLimits(-Inf(nq, 1), Inf(nq, 1));
SL = compile(SL); 

%% Build Actuators
dp.Vb = 225;
dp.Vg = 0;

nact = nu;
sl_actuators = HamrActuators(nact, {'FLsact', 'FLlact'}, [1; 1], dp);

%% Connect system

% connections from actuators to SL
sl_actuators = sl_actuators.setOutputFrame(SL.getInputFrame());
connection1(1).from_output = sl_actuators.getOutputFrame();
connection1(1).to_input = SL.getInputFrame();

% connections from SL to actuators
sl_out = SL.getOutputFrame();
act_in = sl_actuators.getInputFrame();
act_in = act_in.replaceFrameNum(2, sl_out.getFrameByName('ActuatorDeflection'));
sl_actuators = sl_actuators.setInputFrame(act_in);
%
connection2(1).from_output = sl_out.getFrameByName('ActuatorDeflection');
connection2(1).to_input = act_in.getFrameByName('ActuatorDeflection');

% mimo inputs
input_select(1).system = 1;
input_select(1).input = act_in.getFrameByName('DriveVoltage');

% mimo outputs
output_select(1).system = 2;
output_select(1).output = sl_out.getFrameByName('FrontLeftSingleLegPosition');
output_select(2).system = 2;
output_select(2).output = sl_out.getFrameByName('FrontLeftSingleLegVelocity');
output_select(3).system = 1;
output_select(3).output = sl_actuators.getOutputFrame();

SLWact = mimoFeedback(sl_actuators, SL, connection1, connection2, ...
    input_select, output_select);

%% Build and Plot Input

tramp = 500; 
Vscale = (dp.Vb - dp.Vg)/3.5;
V = Vscale*idinput([N, SL.getNumInputs()], 'rgs', band)'; 
V(t < 500) = (t(t<500)/tramp).*V(t<500);
V = V + (dp.Vb - dp.Vg)/2 + dp.Vg;
V(V > dp.Vb) = dp.Vb;
V(V < dp.Vg) = dp.Vg;

figure(1); clf;
subplot(2,1,1);
yyaxis left; hold on; plot(t, V(1,:))
subplot(2,1,2);
yyaxis left;  hold on; plot(t, V(2,:))

freq = [-(Fs/2):(1/N/dt):(Fs/2)];

f0 = 0;
ind0 = find(freq > f0, 1, 'first');
ind1 = find(freq > 2*f1, 1, 'first');
Vfft = abs(fftshift((fft(bsxfun(@minus, V, (dp.Vb - dp.Vg)/2 + dp.Vg), [], 2)/size(V,2))));

figure(2); clf;
subplot(2,1,1); 
yyaxis left; hold on; plot(freq(ind0:ind1)*1e3, 2*Vfft(1,ind0:ind1))
subplot(2,1,2); hold on;
yyaxis left; hold on; plot(freq(ind0:ind1)*1e3, 2*Vfft(2,ind0:ind1))

%% IC and  Simulate

Vtraj = PPTrajectory(zoh(t, V));
Vtraj = Vtraj.setOutputFrame(SLWact.getInputFrame);
SLWact_OL = cascade(Vtraj, SLWact);

q0 = zeros(nq,1);
x0 = [q0; 0*q0];

disp('Simulating...'); 
tic; xtraj = simulate(SLWact_OL, [0 T], x0); tlcp = toc;
x = xtraj.eval(t);

%% Estimated Jacobian

act_inputs = [1, 5];    % index of actuator inputs
sfb_inputs = [3; 7];    % index of SFB inputs

% calculate jacobian (1st order)
u = x(act_inputs, :);
y = x(sfb_inputs, :);
J_sfb_a_hat = u'\y';

% check 
yhat = J_sfb_a_hat*u; 

figure(3); clf;
subplot(2,1,1); hold on;
plot(t, y(1,:)); plot(t, yhat(1, :), '--')
yyaxis left; hold on; plot(t, u(1,:));
subplot(2,1,2); hold on;
plot(t, y(2,:)); plot(t, yhat(2, :), '--')
yyaxis left; hold on; plot(t, u(2,:));

%% Analytic Jacobian 

% v_sfb = J_sfb_a * v_a
J_sfb_a = zeros(nq, numel(sfb_inputs), N);
J_sfb_a(act_inputs, :, :) = repmat(eye(2), 1, 1, N);

for i = 1:N
    q = x(1:nq, i);
    qd = x(nq+(1:nv), i); 
    kinsol = SL.doKinematics(q, qd, struct('compute_gradients', true)); 
    [K, dK] = SL.positionConstraints(q);
    dKUnique = dK(SL.valid_loops, :);
    unactuated_dof = find((1:nq ~= act_inputs(1)) &  (1:nq ~= act_inputs(2)))';
    Jc = -dKUnique(:,unactuated_dof)\dKUnique(:, act_inputs);
    J_sfb_a(unactuated_dof, :, i) = Jc;
end

%% Plot J as function of qswing

figure(4); clf; 
subplot(2, 2, 1); hold on;

xlim([-0.3, 0.3]); 
plot(y(1, :), u(1, :), '.');  
ylabel('Swing Angle')
% plot(u(1, :), J_sfb_a_hat(1, 1)*ones(N, 1), 'LineWidth', 2); 

subplot(2, 2, 2); hold on;
plot(y(1, :), u(2, :), '.');
ylabel('Lift Angle')
% plot(u(1, :), J_sfb_a_hat(1, 2)*ones(N, 1), 'LineWidth', 2);

subplot(2, 2, 3); hold on;
plot(y(2, :), u(1, :), '.');  
ylabel('Lift Angle')
xlabel('Swing Deflection')
% plot(u(1, :), J_sfb_a_hat(2, 1)*ones(N, 1), 'LineWidth', 2);  

subplot(2, 2, 4); hold on;
plot(y(2, :), u(2, :), '.');
xlabel('Lift Deflection')
% plot(u(1, :), J_sfb_a_hat(2, 2)*ones(N, 1), 'LineWidth', 2);  

%% Plot J as function of qswing

figure(5); clf; 
subplot(2, 2, 1); hold on;

title("J_{11}"); xlim([-0.3, 0.3]); 
plot(u(1, :), squeeze(J_sfb_a(sfb_inputs(1), 1, :)), '.');  
plot(u(1, :), J_sfb_a_hat(1, 1)*ones(N, 1), 'LineWidth', 2); 

subplot(2, 2, 2); hold on;
title("J_{12}"); xlim([-0.3, 0.3]); 
plot(u(1, :), squeeze(J_sfb_a(sfb_inputs(1), 2, :)), '.');
plot(u(1, :), J_sfb_a_hat(1, 2)*ones(N, 1), 'LineWidth', 2); 

subplot(2, 2, 3); hold on;
title("J_{21}"); xlim([-0.3, 0.3]); 
plot(u(1, :), squeeze(J_sfb_a(sfb_inputs(2), 1, :)), '.'); 
plot(u(1, :), J_sfb_a_hat(2, 1)*ones(N, 1), 'LineWidth', 2); 


subplot(2, 2, 4); hold on;
title("J_{22}"); xlim([-0.3, 0.3]); 
plot(u(1, :), squeeze(J_sfb_a(sfb_inputs(2), 2, :)), '.'); 
plot(u(1, :), J_sfb_a_hat(2, 2)*ones(N, 1), 'LineWidth', 2);

%% Plot J as function of qlift
figure(6); clf; 
subplot(2, 2, 1); hold on;

title("J_{11}"); xlim([-0.3, 0.3]); 
plot(u(2, :), squeeze(J_sfb_a(sfb_inputs(1), 1, :)), '.');  
plot(u(2, :), J_sfb_a_hat(1, 1)*ones(N, 1), 'LineWidth', 2); 

subplot(2, 2, 2); hold on;
title("J_{12}"); xlim([-0.3, 0.3]); 
plot(u(2, :), squeeze(J_sfb_a(sfb_inputs(1), 2, :)), '.');
plot(u(2, :), J_sfb_a_hat(1, 2)*ones(N, 1), 'LineWidth', 2); 

subplot(2, 2, 3); hold on;
title("J_{21}"); xlim([-0.3, 0.3]); 
plot(u(2, :), squeeze(J_sfb_a(sfb_inputs(2), 1, :)), '.'); 
plot(u(2, :), J_sfb_a_hat(2, 1)*ones(N, 1), 'LineWidth', 2); 

subplot(2, 2, 4); hold on;
title("J_{22}"); xlim([-0.3, 0.3]); 
plot(u(2, :), squeeze(J_sfb_a(sfb_inputs(2), 2, :)), '.'); 
plot(u(2, :), J_sfb_a_hat(2, 2)*ones(N, 1), 'LineWidth', 2);

