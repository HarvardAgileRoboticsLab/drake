clear; clc; close all;
warning('off', 'MATLAB:nargchk:deprecated');

%% Load Rigid Body Model

urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2.urdf');

% options
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.z_inactive_guess_tol = 0.1;
options.use_bullet = false;

% options to change

dt = 0.2;
FREQ = 0.01;
gait = 'TROT';
NCYC = 25;
NPTS = 100;             % number of pts/cycle in desired traj
RAMPCYC = 10;           % number of cycles to ramp
DC = 50;                % duty cycle for swing
DL = 50;
LIFTAMP = 0.15;         % lift actuator motion (mm)
SWINGAMP = 0.175;       % swing actuator motion (mm)
MU = 0.51;                                  % this needs to be manually set in @RigidBodyManipulator/ContactConstraints
TYPE = 1;
SAVE_FLAG = 0;

ISFLOAT = true; % floating (gnd contact) or in air (not floating)

% general RBM options
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

% hamr options
options.stiffness_mult = [1, 1, 1, 1]';     %(swing act, lift act, swing flex, lift flex)

% Build robot + visualizer
hamr = HamrTSRBM(HamrRBM(urdf, options), dt, options);
nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nu = hamr.getNumInputs();
qa = hamr.getActuatedJoints();

v = hamr.constructVisualizer();

%% Build Actuators
dp.Vb = 225;
dp.Vg = 0;

nact = 8;
hr_actuators = HamrActuators(nact, {'FLsact', 'FLlact', 'RLsact', 'RLlact', ...
    'FRsact', 'FRlact', 'RRsact', 'RRlact'}, [], dp);

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

%% Simulate closed loop

params.NCYC = NCYC;
params.RAMPCYC = RAMPCYC;
params.LIFTAMP = LIFTAMP;
params.SWINGAMP = SWINGAMP;
params.GAIT = gait;
params.NPTS = NPTS;
params.TYPE = TYPE;
params.STIFFNESS_MULT = options.stiffness_mult;
params.MU = MU;
params.FREQ = FREQ;
params.DC_SWING = DC;
params.DC_LIFT = DL;
params.X0 = x0;


[tt_sol, xx_sol, vv_sol, err_sol, xtrajd] = HuersiticTrajCLSim(hamrWact, hamr, hr_actuators, params);
ttd = xtrajd.getBreaks();
xxd = xtrajd.eval(ttd);

%% Plotting 

act_dof = hamr.getActuatedJoints();
ndof = 2*hamr.getNumPositions();
title_str = {'Front Left Swing', 'Front Left Lift', ...
    'Rear Left Swing', 'Rear Left Lift', ...
    'Front Right Swing', 'Front Right Lift', ...
    'Rear Rear Swing', 'Rear Rear Lift'};

figure(2); clf; hold on;
for i = 1:numel(act_dof)
    subplot(4,2,i); hold on; title(title_str{i})
%     plot(tt_sol, vv_sol_OL(i,:));
    plot(tt_sol, vv_sol(i,:));
%     legend('OL Inputs', 'CL Inputs');
    %     ylim([dp.Vg, dp.Vb])
end

figure(3); clf; hold on;
for i = 1:numel(act_dof)
    subplot(4,2,i); hold on; title(title_str{i})
%     plot(tt_sol, xx_sol(act_dof(i), :));
    plot(tt_sol, xx_sol(act_dof(i), :));
    plot(ttd, xxd(act_dof(i), :));
    legend('CL Act Defl', 'Desired Act Defl');
    %     ylim([-SWINGAMP, SWINGAMP])
end

figure(4); clf; hold on;
for i = 1:numel(act_dof)
    subplot(4,2,i); hold on; title(title_str{i})
%     plot(tt_sol, xx_sol(nq+act_dof(i), :));
    plot(tt_sol, xx_sol(nq+act_dof(i), :));
    plot(ttd, xxd(nq+act_dof(i), :));
    legend('CL Act Vel', 'Desired Act Vel');
end

fpopt.loc = 'foot';
fpopt.base_frame = 'World';
pfCL = zeros([numel(tt_sol), size(hamr.HamrRBM.FOOT_POS')]);
vfCL = zeros([numel(tt_sol), size(hamr.HamrRBM.FOOT_POS')]);
nc = numel(hamr.HamrRBM.LEG_NAME);

for i = 1:numel(tt_sol)
    [pfCL(i,:,:), Ji] = hamr.HamrRBM.getFootPosition(xx_sol(1:nq,i), xx_sol(nq+(1:nv),i), fpopt);
    for j = 1:numel(Ji)
        vfCL(i, :,j) = Ji{j}*xx_sol(nq+(1:nv),i);
    end
end

% figure(5); clf; hold on;
% for i = 1:nc
%     subplot(nc/2,2,i); hold on;
%     xlabel('Foot X'); ylabel('Foot VX'); title(hamr.HamrRBM.LEG_NAME{i})
%     plot(pfCL(:,1,i) - pfCL(1,1,i) , vfCL(:,1,i))
%     axis equal;
% end
%
% figure(6); clf; hold on;
% for i = 1:nc
%     subplot(nc/2,2,i); hold on;
%     xlabel('Foot Z'); ylabel('Foot VZ'); title(hamr.HamrRBM.LEG_NAME{i})
%     plot(pfCL(:,3,i), vfCL(:,3,i))
%     axis equal;
% end

figure(7); clf; hold on;
for i = 1:nc
    subplot(nc/2,2, i); hold on;
    xlabel('Foot X'); ylabel('Foot Z'); title(hamr.HamrRBM.LEG_NAME{i})
    plot(pfCL(:,1,i) - pfCL(1,1,i), pfCL(:,3,i)- pfCL(1,3,i));
    axis equal;
    %     xlim([-2, 2])
    %     ylim([-2, 2])
end

if options.floating
    figure(8); clf;
    title_str = {'Com-X(mm)', 'Com-Y(mm)', 'Com-Z(mm)', 'Roll(deg)', 'Pitch(deg)', 'Yaw(deg)'};
    for i = 1:6
        subplot(2,3,i); hold on; ylabel(title_str(i), 'FontSize', 18)
        xlabel('Time(s)', 'FontSize', 18)
        if i > 3
%             plot(tt_sol*1e-3, rad2deg(xx_sol(i,:)), 'LineWidth', 1.5);
            plot(tt_sol*1e-3, rad2deg(xx_sol(i,:)), 'LineWidth', 1.5);
        else
            plotyy(tt_sol*1e-3, xx_sol(i,:),tt_sol*1e-3, xx_sol(nq+i,:));
        end
    end
end

tilefigs;

%% Playback
xtraj_scaled = DTTrajectory(tt_sol*1e-3, xx_sol(1:nq,:));
xtraj_scaled = xtraj_scaled.setOutputFrame(v.getInputFrame());
options.slider = true;
v.playback(xtraj_scaled, options);

