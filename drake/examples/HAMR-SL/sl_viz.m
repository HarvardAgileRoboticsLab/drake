clear; clc; close all;
format long

%% User inputs

fd = 0.1;         % drive frequency (Hz)
tsim = 50/fd; 

% ramp
tramp = 10/fd;

dp.Vb = 50;
dp.Vg = 0;

%% Load Rigid Body

pkg_path = '/home/nddoshi/dev/drake/drake/examples/HAMR-SL/';
urdf = [pkg_path, 'SL_assem6/urdf/SL_scaled.urdf'];

% options
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.z_inactive_guess_tol = 1;
options.dt = 1/(50*fd);  %time step
options.use_bullet = false;  
options.floating = false;
options.collision = false;
options.terrain = []; 

x0 = zeros(16,1); 

% Build robot + visualizer
hamrsl = HamrSL(urdf, options);
hamrsl = compile(hamrsl);
v = hamrsl.constructVisualizer();


%% Build Actuators

nact = 2;
hr_actuators = HamrActuators(nact, {'FLsact', 'FLlact'}, [1; 1;], dp);

%% Connect system

% connections from actuators to hamr
hr_actuators = hr_actuators.setOutputFrame(hamrsl.getInputFrame());
connection1(1).from_output = hr_actuators.getOutputFrame();
connection1(1).to_input = hamrsl.getInputFrame();

% connections from hamr to actuators
hamr_out = hamrsl.getOutputFrame();
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

hamrWact = mimoFeedback(hr_actuators, hamrsl, connection1, connection2, ...
    input_select, output_select);

%% Build (open-loop) control input

t = 0:options.dt:tsim;

%TROT
Vact = [0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t + pi/2);            % FLswing
    0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t)];                      % FLlift
ramp = t/tramp; ramp(t >= tramp) = 1;

Vact = bsxfun(@times, ramp, Vact) + 0.5*(dp.Vb - dp.Vg);
u = PPTrajectory(foh(t, Vact));
u = setOutputFrame(u, hamrWact.getInputFrame());

figure(1); clf;
plot(t, Vact(1,:), t, Vact(2,:), '--');
legend('Swing Drive', 'Lift Drive')

%% Simulate Open loop

hamr_OL = cascade(u, hamrWact);
x0_hat = hamrsl.getManipulator.positionConstraints(x0); 
[tf, err_str] = valuecheck(positionConstraints(hamrsl,x0),zeros(72,1),1e-8);
tf = 1; 
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

%% Plotting (last 10 cycles)

tt = xtraj_scaled.getBreaks();
yy = xtraj_scaled.eval(tt);

ind_s = find(tt > tsim*1e-3 - 

act_dof = hamrsl.getActuatedJoints();
ndof = hamrsl.getNumDiscStates();
title_str = {'Front Left Swing', 'Front Left Lift', ...
    'Rear Left Swing', 'Rear Left Lift', ...
    'Front Right Swing', 'Front Right Lift', ...
    'Rear Rear Swing', 'Rear Rear Lift'};

figure(2); clf; hold on;
for i = 1:numel(act_dof)
    subplot(2,1,i); hold on; title(title_str(i))
    yyaxis left; hold on; plot(tt, yy(ndof+i,:), 'b')
%     yyaxis left; plot(tt, yy(act_dof(i), :)*1e3); 
%     yyaxis right; plot(tt, Vact(i,:)); 
%     legend('Deflection', 'Force')
    yyaxis right; hold on; plot(tt, yy(act_dof(i), :)*1e3, 'r--', ...
        t/1000, Vact(i,:) - mean(Vact(i,:)), 'r')
    legend('Force', 'Deflection', 'Drive')
end

lp_b = [0, 7.540, -11.350];

lp_g = zeros([numel(t), size(lp_b')]); 

legs = {'FLL4'}; 
% w = warning('off','Drake:RigidBodyManipulator:WeldedLinkInd');
for j = 1:numel(tt)
    q = yy(1:ndof/2, j);
    qd = yy(ndof/2+1: ndof, j);
    kinsol = hamrsl.doKinematics(q, qd);
    for i = 1:size(lp_b,1)
        lp_g(j,:,i) = hamrsl.forwardKin(kinsol, hamrsl.findLinkId(legs{i}), lp_b(i,:)'); 
    end
end
% warning(w);

figure(3); clf; hold on;
for i = 1:size(lp_b,1)
    title(legs{i});
    plot((lp_g(:,1,i) - mean(lp_g(:,1,i))), ...
        (lp_g(:,3,i) - mean(lp_g(:,3,i))))    
    axis equal; 
end

