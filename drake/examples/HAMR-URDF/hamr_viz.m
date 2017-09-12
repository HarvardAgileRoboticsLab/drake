clear; clc; close all;
% format long

%% Load Rigid Body

urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaled.urdf');

% options
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.z_inactive_guess_tol = 1;
options.dt = 0.1;  %time step
options.use_bullet = false;  

% floating (gnd contact) or in air (not floating)
ISFLOAT = false; 
% ISFLOAT = true; 
if ISFLOAT
    options.floating = ISFLOAT;
    options.collision = ISFLOAT;
    x0 = zeros(76,1); x0(3) = 12;
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
v = hamr.constructVisualizer();

%% Build Actuators
dp.Vb = 95;
dp.Vg = 0;
% 
nact = 8;
hr_actuators = HamrActuators(nact, {'FLsact', 'FLlact', 'RLsact', 'RLlact', ...
    'FRsact', 'FRlact', 'RRsact', 'RRlact'}, [1; 1; -1; -1; 1; 1; -1; -1], dp);
% % hr_actuators = HamrActuators(nact, {'FLsact', 'FLlact', 'RLsact', 'RLlact'}, dp);
% % hr_actuators = HamrActuators(nact, {'FLsact', 'FLlact'}, dp);

%% Connect system

% connections from actuators to hamr
hr_actuators = hr_actuators.setOutputFrame(hamr.getInputFrame());
connection1(1).from_output = hr_actuators.getOutputFrame();
connection1(1).to_input = hamr.getInputFrame();

% connections from hamr to actuators
hamr_out = hamr.getOutputFrame();
act_in = hr_actuators.getInputFrame();
act_in = act_in.replaceFrameNum(2, hamr_out.getFrameByName('ActuatorDeflectionandRate'));
hr_actuators = hr_actuators.setInputFrame(act_in);
% 
connection2(1).from_output = hamr_out.getFrameByName('ActuatorDeflectionandRate');
connection2(1).to_input = act_in.getFrameByName('ActuatorDeflectionandRate');

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

fd = 0.065;         % drive frequency (Hz)
tsim = 0.2e3; 

t = 0:options.dt:tsim;

% Vact = [0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t + pi/2);            % FLswing
%     0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t - deg2rad(60));                       % FLlift
%     0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t + 3*pi/2);                % RLSwing
%     0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t - deg2rad(60));                  % RLLift
%     0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t + pi/2 + deg2rad(120));                % FRswing
%     0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t + deg2rad(60));                       % FRlift
%     0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t + 3*pi/2 - deg2rad(120));              % RRSwing
%     0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t + deg2rad(60))];                      % RRLift

%TROT
Vact = [0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t + pi/2);            % FLswing
    0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t);                       % FLlift
    0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t + pi/2);              % RLSwing
    0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t - pi);                       % RLLift
    0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t + pi/2);                % FRswing
    0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t);                       % FRlift
    0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t + pi/2);              % RRSwing
    0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t - pi)];                      % RRLift

% ramp
tramp = 3/fd;
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
    v.playback(xtraj_scaled, options);
else
    disp('invalid initial condition...')
end

%% Plotting
tt = xtraj.getBreaks();
yy = xtraj.eval(tt);

act_dof = hamr.getActuatedJoints();
ndof = hamr.getNumDiscStates();
title_str = {'Front Left Swing', 'Front Left Lift', ...
    'Rear Left Swing', 'Rear Left Lift', ...
    'Front Right Swing', 'Front Right Lift', ...
    'Rear Rear Swing', 'Rear Rear Lift'};

figure(2); clf; hold on;
for i = 1:numel(act_dof)
    subplot(4,2,i); hold on; title(title_str(i))
    yyaxis left; hold on; plot(tt, yy(ndof+i,:), 'b')
%     yyaxis left; plot(tt, yy(act_dof(i), :)*1e3); 
%     yyaxis right; plot(tt, Vact(i,:)); 
%     legend('Deflection', 'Force')
    yyaxis right; hold on; plot(tt, yy(act_dof(i), :)*1e3, 'r--', ...
        t/1000, Vact(i,:) - mean(Vact(i,:)), 'r')
    legend('Force', 'Deflection', 'Drive')
end

lp_b = [0, 7.540, -11.350;
    0, 7.540, -11.350;
    0, -7.540, -11.350;
    0, -7.540, -11.350]; 

lp_g = zeros([numel(t), size(lp_b')]); 

legs = {'FLL4', 'RLL4', 'FRL4', 'RRL4'}; 

for j = 1:numel(tt)
    q = yy(1:ndof/2, j);
    qd = yy(ndof/2+1: ndof, j);
    kinsol = hamr.doKinematics(q, qd);
    for i = 1:size(lp_b,1)
        lp_g(j,:,i) = hamr.forwardKin(kinsol, hamr.findLinkId(legs{i}), lp_b(i,:)'); 
    end
end

figure(3); clf; hold on;
for i = 1:size(lp_b,1)
%     subplot(2,2,i); hold on; title(legs{i});
    plot((lp_g(:,1,i) - mean(lp_g(:,1,i))), ...
        (lp_g(:,3,i) - mean(lp_g(:,3,i))))  
    axis equal; 
%     axis([-2.5, 2.5, -2.5, 2.5])
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
    
    figure(5); clf; hold on;
    plot(tt, phi);
end
