clear; clc; close all;
% global kl_traj jl_traj c_traj beta_traj psi_traj eta_traj
save_dir = '~/Dropbox/CurrentWork/FrictionTrajOpt/MatFiles/forXPC/';

%% Load Rigid Body

urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2.urdf');

% kl_traj = [];
% jl_traj = [];
% c_traj = [];
% beta_traj = [];
% psi_traj = [];
% eta_traj = [];

% options
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.z_inactive_guess_tol = 0.1;
options.use_bullet = false;

% options to change
options.dt = 1;
% options.mu = 0.6;
gait = 'TROT';
SAVE_FLAG = 1;
ISFLOAT = false; % floating (gnd contact) or in air (not floating)

if ISFLOAT
    options.floating = ISFLOAT;
    options.collision = ISFLOAT;
    x0 = zeros(76, 1); x0(3) = 12.69;
    options.terrain = RigidBodyFlatTerrain();
    
else
    options.floating = ISFLOAT;
    options.collision = ISFLOAT;
    x0 = zeros(64, 1);
    options.terrain = [];
end

% Build robot + visualizer
hamr = HamrTSRBM(urdf, options);
hamr = compile(hamr);
nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
% v = hamr.constructVisualizer();
% v.inspector(x0);

%% Build Actuators
dp.Vb = 150;
dp.Vg = 0;

nact = 8;
hr_actuators = HamrActuators(nact, {'FLsact', 'FLlact', 'RLsact', 'RLlact', ...
    'FRsact', 'FRlact', 'RRsact', 'RRlact'}, [-1; -1; 1; 1; -1; -1; 1; 1], dp);

%% Connect system

% connections from actuators to hamr
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
output_select(3).system = 1;
output_select(3).output = hr_actuators.getOutputFrame();

hamrWact = mimoFeedback(hr_actuators, hamr, connection1, connection2, ...
    input_select, output_select);

%% Build (open-loop) control input

fd = 0.01;         % drive frequency (Hz)
tsim = 1e3;

t = 0:options.dt:tsim;

switch gait
    case 'TROT'
        Vact = [0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t + pi/2);            % FLswing
            0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t);                       % FLlift
            0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t + 3*pi/2);              % RLSwing
            0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t);                       % RLLift
            0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t + pi/2);                % FRswing
            0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t);                       % FRlift
            0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t + 3*pi/2);              % RRSwing
            0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t)];                      % RRLift
    case 'PRONK'
        Vact = [0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t + pi/2);            % FLswing
            0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t);                       % FLlift
            0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t + pi/2);                % RLSwing
            0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t + pi);                  % RLLift
            0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t - pi/2);                % FRswing
            0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t + pi);                       % FRlift
            0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t - pi/2);                % RRSwing
            0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t)];                 % RRLift
    otherwise
        Vact = zeros(8, numel(t));
end

% ramp
tramp = 1/fd;
ramp = t/tramp; ramp(t >= tramp) = 1;

Vact = bsxfun(@times, ramp, Vact) + 0.5*(dp.Vb - dp.Vg);
u = PPTrajectory(foh(t, Vact));
u = setOutputFrame(u, hamrWact.getInputFrame());

figure(1); clf;
plot(t, Vact(1,:), t, Vact(2,:), '--');
legend('Swing Drive', 'Lift Drive')

%% Simulate Open loop

hamr_OL = cascade(u, hamrWact);
nQ = hamr.getNumPositions();
x0_hat = hamr.getManipulator.positionConstraints(x0);
[tf, err_str] = valuecheck(positionConstraints(hamr,x0),zeros(72,1),1e-6);

% x0_hat = hamr.positionConstraints(x0(1:nQ));
% [tf, err_str] = valuecheck(positionConstraints(hamr,x0(1:nQ)),zeros(72,1),1e-6);

% tf = 1;
if tf
    disp('Valid initial condition: simulating...')
    tic;
    xtraj = simulate(hamr_OL, [0 tsim], x0);
    tlcp = toc;
    xtraj_scaled = DTTrajectory(xtraj.getBreaks()*1e-3, xtraj.eval(xtraj.getBreaks()));
    xtraj_scaled = xtraj_scaled.setOutputFrame(xtraj.getOutputFrame());
    fprintf('It took %fs to simulate %fs of realtime. \nThats %fx \n', ...
        tlcp, tsim/1000, 1000*tlcp/tsim)
    options.slider = true;
    %     xtraj.tt = xtraj.tt/1000;
%     v.playback(xtraj_scaled, options);
else
    disp('invalid initial condition...')
end

%% Plotting
tt = xtraj.getBreaks();
yy = xtraj.eval(tt);
xx = yy(1:2*nQ,:);
uu = yy(2*nQ+1:end,:);

act_dof = hamr.getActuatedJoints();
swing_dof = act_dof(1:2:end); 
title_str = {'Front Left Swing', 'Front Left Lift', ...
    'Rear Left Swing', 'Rear Left Lift', ...
    'Front Right Swing', 'Front Right Lift', ...
    'Rear Rear Swing', 'Rear Rear Lift'};

figure(2); clf; hold on;
for i = 1:numel(act_dof)
    subplot(4,2,i); hold on; title(title_str(i))
    yyaxis left; hold on; plot(tt, uu(i,:))
    %     yyaxis left; plot(tt, yy(act_dof(i), :)*1e3);
    %     yyaxis right; plot(tt, Vact(i,:));
    %     legend('Deflection', 'Force')
    yyaxis right; hold on; plot(tt, yy(act_dof(i), :)*1e3)
    legend('Force', 'Deflection')
end

lp_b = [0, 8.58, -6.00;
    0, 8.58, -6.00;
    0, -8.58, -6.00;
    0, -8.58, -6.00];

pf0 = 0*lp_b;
pfB = zeros([numel(t), size(lp_b')]);
vfB = zeros([numel(t), size(lp_b')]);

legs = {'FLL4', 'RLL4', 'FRL4', 'RRL4'};

kinsol0 = hamr.doKinematics(x0(1:nq), x0(nq+(1:nv))); 
for i = 1:size(lp_b,1)
    pf0(i,:) = hamr.forwardKin(kinsol0, hamr.findLinkId(legs{i}), lp_b(i,:)')';
end

for i = 1:size(lp_b,1)
    for j = 1:numel(tt)
        q = yy(1:nq, j);
        qd = yy(nq+(1:nv), j);
        kinsol = hamr.doKinematics(q, qd, struct('compute_gradients', true));
        [pfB_ji, dpf_dq] = hamr.forwardKin(kinsol, hamr.findLinkId(legs{i}), lp_b(i,:)');
        pfB(j,:,i) = pfB_ji - pf0(i,:)';
        vfB(j,:,i) = dpf_dq*qd; 
    end
end

figure(3); clf; hold on;
for i = 1:size(lp_b,1)
%         subplot(2,2,i); hold on; title(legs{i});
    plot(pfB(:,1,i), pfB(:,3,i)); 
    axis equal;
%         axis([-2.5, 2.5, -2.5, 2.5])
end
legend(legs)

if ISFLOAT
    figure(4); clf;
    title_str = {'com x', 'com y', 'com z', 'roll', 'pitch', 'yaw'};
    for i = 1:6
        subplot(3,2,i); hold on; title(title_str(i))
        yyaxis left; hold on; plot(tt, yy(i,:))
        yyaxis right; hold on; plot(tt, yy(i+ndof/2, :))
    end
    
    contact_opt.use_bullet = false;
    phi = zeros(4, numel(tt));
    for i = 1:numel(tt)
        q = yy(1:ndof/2, i);
        qd = yy(ndof/2+1:ndof, i);
        kinsol = hamr.doKinematics(q, qd);
        phi(:,i) = hamr.contactConstraints(kinsol, false, contact_opt);
    end
    
    figure(5); clf;
    for i = 1:size(phi, 1)
        subplot(2,2,i); hold on;
        yyaxis right; plot(tt, phi(i,:)); ylabel('Leg Height')
        yyaxis left; plot(tt, c_traj(i,:)); ylabel('Force')
        plot(tt, beta_traj(4*(i-1)+(1:4), :), 'k');
    end
end

%% saving

nc = tsim*fd; 
nCF0 = find(tt >= (nc-1)/fd, 1, 'first');
nCF1 = numel(tt); 
% 
% % zero out starting x and y pos
% xx_sim = sim_traj.xx(:, nCF0:nCF1);
% xx_sim0 = 0*xx_sim(:,1); xx_sim0(1:2) = xx_sim(1:2, 1); 
% 
% xx_sim = bsxfun(@minus, xx_sim, xx_sim0); 
tt = (tt(nCF0:nCF1) - tt(nCF0))/1e3; 
xx = xx(:, nCF0:nCF1); 
uu = Vact(:, nCF0:nCF1); 
pfB = pfB(nCF0:nCF1, :, :); 
vfB = vfB(nCF0:nCF1, :, :); 

if SAVE_FLAG
    fname = ['NC_', gait, '_', num2str(dp.Vb) 'N_', num2str(1e3*fd), 'Hz' '.mat'];
    disp(['Saving as ', fname]);
    save([save_dir, fname], 'tt', 'xx', 'uu', 'pfB', 'vfB');
end
