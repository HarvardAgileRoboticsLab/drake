clear; clc; close all;
save_dir = '~/Dropbox/GaitRecoveryandTransition/LinearModelFit/';
%global utraj kltraj 

%utraj = [];
%kltraj = []; 

% Sample Rate and Time
dt = 0.4;
Fs = 1/dt;
T = 2000;
t = 0:dt:T;
N = numel(t);

% White Noise Freq Characteristics
f1 = 70e-3;
band = [0, 2*f1*dt];

%% Build robot

name = 'HAMR_scaledV2';
urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', [name, '.urdf']);

options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = false;
options.collision = false;
options.dt = dt; %0.1; 

% Build robot + visualizer
hamr = HamrTSRBM(HamrRBM(urdf, options), options.dt, options);
nq = hamr.getNumPositions(); 
nv = hamr.getNumVelocities(); 
nu = hamr.getNumInputs(); 
v = hamr.constructVisualizer(); 

% Relax joint limits
hamr = hamr.setJointLimits(-Inf(nq, 1), Inf(nq, 1));
hamr = compile(hamr); 

%% Build Actuators

dp.Vb = 225;
dp.Vg = 0;

nact = 8;
act_names = {'FLsact', 'FLlact', 'RLsact', 'RLlact', ...
    'FRsact', 'FRlact', 'RRsact', 'RRlact'};
hr_actuators = HamrActuators(nact, act_names, [], dp);

% for i = 2:2:nact
%     hr_actuators.dummy_bender(i) = hr_actuators.dummy_bender(i).setCFThickness(0.115); 
% end

%% Connect system

%connections from actuators to hamr
hr_actuators = hr_actuators.setOutputFrame(hamr.getInputFrame());
connection1(1).from_output = hr_actuators.getOutputFrame();
connection1(1).to_input = hamr.getInputFrame();

% connections from hamr to actuators
hamr_out = hamr.getOutputFrame();
act_in = hr_actuators.getInputFrame();
act_in = act_in.replaceFrameNum(2, hamr_out.getFrameByName('ActuatorDeflection'));
hr_actuators = hr_actuators.setInputFrame(act_in);
% 
connection2(1).from_output = hamr_out.getFrameByName('ActuatorDeflection');
connection2(1).to_input = act_in.getFrameByName('ActuatorDeflection');

% mimo inputs
input_select(1).system = 1;
input_select(1).input = act_in.getFrameByName('DriveVoltage');

% mimo outputs
output_select(1).system = 2;
output_select(1).output = hamr_out.getFrameByName('HamrPosition');
output_select(2).system = 2;
output_select(2).output = hamr_out.getFrameByName('HamrVelocity');

hamrWact = mimoFeedback(hr_actuators, hamr, connection1, connection2, ...
    input_select, output_select);

%% Build and Plot Input

tramp = 500; 
Vscale = (dp.Vb - dp.Vg)/2;
V = Vscale*idinput([N, hamr.getNumInputs()], 'rgs', band)'; 
V(t < 500) = (t(t<500)/tramp).*V(t<500);
V = V + (dp.Vb - dp.Vg)/2 + dp.Vg;

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
Vtraj = Vtraj.setOutputFrame(hamrWact.getInputFrame);
hamrWact_OL = cascade(Vtraj, hamrWact);

q0 = zeros(nq,1);
x0 = [q0; 0*q0];

disp('Simulating...'); 
tic; xtraj = simulate(hamrWact_OL, [0 T], x0); tlcp = toc;
x = xtraj.eval(t);

%%  Plot Response
act_dof = hamr.getActuatedJoints(); 
y = x([act_dof; act_dof+nq], :);
yfft = abs(fftshift((fft(y, [], 2)/size(y,2))));

% figure(1);
% subplot(2,1,1); 
% yyaxis right; hold on; plot(t, rad2deg(y(1,:)));
% subplot(2,1,2); hold on;
% yyaxis right; hold on; plot(t, rad2deg(y(2,:)));
% 
% figure(2);
% subplot(2,1,1);  hold on;
% yyaxis right; hold on;plot(freq(ind0:ind1)*1e3, 2*yfft(1,ind0:ind1))
% subplot(2,1,2); hold on;
% yyaxis right; hold on;plot(freq(ind0:ind1)*1e3, 2*yfft(2,ind0:ind1))

%% Convert Voltage to Torque about Leg

% tau = zeros(nu, N);
% 
% actuated_dof = hamr.getActuatedJoints();
% for i = 1:N
%     q = x(1:nq, i);
%     qd = x(nq+(1:nv), i); 
%     kinsol = hamr.doKinematics(q, qd, struct('compute_gradients', true));
% %     [xf, Jf] = SL.forwardKin(kinsol, SL.findLinkId('FLL4'), zeros(3,1), struct('rotation_type', 1));
% %     [H, C, B] = SL.manipulatorDynamics(q, qd); 
% %     F = pinv(Jf')*B*x(nq+nv+(1:nu),i);
% %     tau(:,i) = F(4:end); 
%     [K, dK] = hamr.positionConstraints(q);
%     dKUnique = dK(hamr.valid_loops, :);
%     unactuated_dof = find((1:nq ~= actuated_dof(1)) &  (1:nq ~= actuated_dof(2)))';
%     Jc = -dKUnique(:,unactuated_dof)\dKUnique(:, actuated_dof);
%     Jc_phi = Jc([2, 5], :);    
%     tau(:,i) = Jc_phi\x(nq+nv+(1:nu),i);    
% end
% 
% figure(3); clf;
% subplot(2,1,1); hold on; 
% title('Swing Input')
% plot(t, x(nq+nv+1,:)); 
% plot(t, tau(1,:)); 
% legend('U', 'Tau')
% 
% subplot(2,1,2); hold on; 
% title('Lift Input')
% plot(t, x(nq+nv+2,:)); 
% plot(t, tau(2,:)); 
% legend('U', 'Tau')

disp('Saving...')
save([save_dir, 'sysid_traj_flipped_', name, '_', num2str(f1*1e3), 'Hz_', ...
    num2str(dp.Vb), 'V'], 't', 'V', 'y'); 
