clear; clc; close all;
addpath('../')

urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', ...
    'urdf', 'HAMR_scaledV2.urdf');

% options
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.z_inactive_guess_tol = 0.1;
options.use_bullet = false;

% options to change
options.dt = 0.25;
ISFLOAT = true; % floating (gnd contact) or in air (not floating)
options.floating = ISFLOAT;
options.collision = ISFLOAT;
options.terrain = RigidBodyFlatTerrain();


% Build robot + visualizer
hamr = HamrTSRBM(HamrRBM(urdf,options), options.dt, options);
hamr = compile(hamr);
v = hamr.constructVisualizer();

% state/input dimenisons
nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nx = nq+nv;
nu = hamr.getNumInputs();

%% Build Actuators
dp.Vb = 225;
dp.Vg = 0;
%
nact = 8;
hr_actuators = HamrActuators(nact, {'FLsact', 'FLlact', 'RLsact', 'RLlact', ...
    'FRsact', 'FRlact', 'RRsact', 'RRlact'}, [-1; -1; 1; 1; -1; -1; 1; 1], dp);

%% Define system connections

%connections from actuators to hamr
hr_actuators = hr_actuators.setOutputFrame(hamr.getInputFrame());
connection1(1).from_output = hr_actuators.getOutputFrame();
connection1(1).to_input = hamr.getInputFrame();

% connections from hamr to actuators
hamr_out = hamr.getOutputFrame();
act_in = hr_actuators.getInputFrame();
act_in = act_in.replaceFrameNum(2, hamr_out.getFrameByName('ActuatorDeflection'));
hr_actuators = hr_actuators.setInputFrame(act_in);

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


%% Build hamr with actuators

hamrWact = mimoFeedback(hr_actuators, hamr, connection1, connection2, ...
    input_select, output_select);

%% Build web request system

HWebRequest = HamrWebRequest(hamrWact, hamr, 'true'); 

% mimo outputs
output_select(1).system = 1;
output_select(1).output = hamrWact.getOutputFrame.getFrameByName('HamrPosition');
output_select(2).system = 1;
output_select(2).output = hamrWact.getOutputFrame.getFrameByName('HamrVelocity');
output_select(3).system = 2;
output_select(3).output = HWebRequest.getOutputFrame();

hamrWact_CL = mimoFeedback(hamrWact, HWebRequest, [], [], [], output_select); 

%% Simulate
T = 150;
tt = 0:options.dt:T;

x0 = zeros(76, 1);
x0(3) = 12.69;
% x0(1:6) = [-6.12311548e-02; 1.36816396e-01; 1.24401108e+01; 
%     -2.80173663e-03; -4.19730058e-03; -3.10018578e-03];

[ytraj, xtraj] = simulate(hamrWact_CL, [0, T], x0);

xx = xtraj.eval(tt); 
yy = ytraj.eval(tt); 

% save('ppohamrfull4', 'tt', 'xx', 'yy')

%% comparison 
% 
% hfull_data = load('ppohamrfull.mat');


xhip = xx(HWebRequest.ihip, 1:end-1);

q = [xx(1:6, 1:end-1); xhip(HWebRequest.perm, :)];
qdot = [xx(HWebRequest.nq + (1:6), 1:end-1); ...
    xhip(HWebRequest.nu + HWebRequest.perm, :) ];
V = yy(HWebRequest.nq + HWebRequest.nv + HWebRequest.perm, :);

save('ppotraj_new_5_hfull', 'q', 'qdot', 'V');


% ppo_dat = load('ppotraj3.mat');
% qppo = ppo_dat.q(:, 1:14)'; 
% qdot_ppo = ppo_dat.q(:, 15:end)'; 
% V_ppo = ppo_dat.u';


figure(1); clf;
for i = 1:14
    subplot(7, 2, i); hold on;
    title(num2str(i))
    plot(q(i, :))
%              plot(qppo(i, :)); 
%     legend('
end

figure(2); clf;
for i = 1:14
    subplot(7, 2, i); hold on;
    title(num2str(i))
        plot(qdot(i, :))
%     plot(qdot_ppo(i, :)); 
%     legend('
end

figure(3); clf;
for i = 1:8
    subplot(4, 2, i); hold on;
    title(num2str(i))
    plot(V(i, :))
%     plot(V_ppo(i, :), '--'); 
%     legend('
end
% xx = xtraj_sim_CL.eval(tt);

%% Playback
xtraj_scaled = DTTrajectory(tt*1e-3, xx);
xtraj_scaled = xtraj_scaled.setOutputFrame(xtraj.getOutputFrame());
options.slider = true;
v.playback(xtraj_scaled, options);

