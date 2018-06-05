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
options.dt = dt;
gait = 'TROT';
DC = 50;                % duty cycle for swing
DL = 50; 
LIFTAMP = 0.15;         % lift actuator motion (mm)     
SWINGAMP = 0.175;       % swing actuator motion (mm)
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


%% Build (open-loop) control input

Vff = 0; %100;          % feed forward voltage
fd = 0.01;              % drive frequency (kHz)
NCYC = 3;
tsim = NCYC/fd;

t = 0:options.dt:tsim;

switch gait    
    case 'TROT'
        vv = [(Vff/2)*sin(2*pi*fd*t);            % FLswing
            (Vff/2)*sin(2*pi*fd*t - pi/2 );                     % FLlift
            (Vff/2)*sin(2*pi*fd*t - pi);              % RLSwing
            (Vff/2)*sin(2*pi*fd*t - pi/2);                     % RLLift
            (Vff/2)*sin(2*pi*fd*t);              % FRswing
            (Vff/2)*sin(2*pi*fd*t - pi/2);                     % FRlift
            (Vff/2)*sin(2*pi*fd*t);              % RRSwing
            (Vff/2)*sin(2*pi*fd*t + pi/2)];                    % RRLift
    case 'PRONK'
        vv = [(Vff/2)*sin(2*pi*fd*t + pi/2);            % FLswing
            (Vff/2)*sin(2*pi*fd*t);                     % FLlift
            (Vff/2)*sin(2*pi*fd*t - pi/2);              % RLSwing
            (Vff/2)*sin(2*pi*fd*t);                     % RLLift
            (Vff/2)*sin(2*pi*fd*t - pi/2);              % FRswing
            (Vff/2)*sin(2*pi*fd*t + pi);                % FRlift
            (Vff/2)*sin(2*pi*fd*t + pi/2);              % RRSwing
            (Vff/2)*sin(2*pi*fd*t + pi)];               % RRLift
    otherwise
        vv = zeros(8, numel(t));
end

% ramp
tramp = 1/fd;
ramp = t/tramp; ramp(t >= tramp) = 1;

vv = bsxfun(@times, ramp, vv) + dp.Vb/2;
vtraj = PPTrajectory(foh(t, vv));
vtraj = setOutputFrame(vtraj, hamrWact.getInputFrame());

%% Generate Hueristic Trajectories in body frame

% generate trajectory 
NPTS = 100;
[traj, brkVal] = GenerateHueristicActTraj(gait, fd, DC, DL, NPTS, TYPE); 
tcyc = traj(:, 1); 
qcyc = traj(:,2:end);

% rescale
qcyc(:, 1:2:end) = SWINGAMP*qcyc(:, 1:2:end);
qcyc(:, 2:2:end) = LIFTAMP*qcyc(:, 2:2:end);
 
% Finite diff to find v_cyc
hhd = (1/fd)/NPTS; 
e = 1./(2*hhd*ones(NPTS,1));
D1 = spdiags([-e, 0*e, e], [-1, 0, 1], NPTS, NPTS);
D1(1,1) = -3/2/hhd; D1(1,2) = 2/hhd; D1(1,3) =-1/2/hhd;
D1(NPTS,NPTS) = 3/2/hhd; D1(NPTS,NPTS-1) = -2/hhd; D1(NPTS,NPTS-2) = 1/2/hhd;

vcyc = D1*qcyc;

tlabel = {'FL', 'RL', 'FR', 'RR'};
figure(1); clf;
for j =1:nact/2
    subplot(nact/4, 2, j); hold on; 
    plot(tcyc, qcyc(:,2*j-1)); 
    plot(tcyc, qcyc(:,2*j));
    legend('Swing', 'Lift')
end

% create trajectory 
ttd = 0:hhd:(tsim-hhd); 
xxd = zeros(nq+nv, NPTS*NCYC); 
xxa = bsxfun(@times, interp1(t, ramp, ttd), repmat([qcyc, vcyc], NCYC, 1)'); 
xxd([qa; nq+qa], :) = xxa;
xtrajd = PPTrajectory(foh(ttd, xxd)); 

%% Simulate Open loop

hamr_OL = cascade(vtraj, hamrWact);
xtraj_sim = simulate(hamr_OL, [0 tsim], x0);

%% Simulate Closed Loop

% sets all errors to ~ 1
Qpos = 2*(2/(LIFTAMP + SWINGAMP))^2;
Qvel = 0; %0.005*(2/(2*pi*fd)/(LIFTAMP + SWINGAMP))^2;
rho = 1*(1/dp.Vb)^2;

tracking_opt.ctype = 'actlqr';
tracking_opt.rho = rho;
tracking_opt.Qpos = Qpos;
tracking_opt.Qvel = Qvel;

LegTracking = HAMRLegTracking(hamrWact, hamr, hr_actuators, vtraj, xtrajd, tracking_opt);

% mimo outputs
output_select(1).system = 1;
output_select(1).output = hamrWact.getOutputFrame.getFrameByName('HamrPosition');
output_select(2).system = 1;
output_select(2).output = hamrWact.getOutputFrame.getFrameByName('HamrVelocity');
output_select(3).system = 2;
output_select(3).output = LegTracking.getOutputFrame();

hamr_CL = mimoFeedback(hamrWact, LegTracking, [], [], [], output_select);
xtraj_sim_CL = simulate(hamr_CL, [0, tsim], x0);

%% Plotting

tt_sol = xtraj_sim.getBreaks();
yy_sol_OL = xtraj_sim.eval(tt_sol);
yy_sol_CL = xtraj_sim_CL.eval(tt_sol);

xx_sol_OL = yy_sol_OL(1:nq+nv, :);
xx_sol_CL = yy_sol_CL(1:nq+nv, :);

vv_sol_OL = vtraj.eval(tt_sol);
vv_sol_CL = yy_sol_CL(nq+nv+(1:nu),:);

act_dof = hamr.getActuatedJoints();
ndof = 2*hamr.getNumPositions();
title_str = {'Front Left Swing', 'Front Left Lift', ...
    'Rear Left Swing', 'Rear Left Lift', ...
    'Front Right Swing', 'Front Right Lift', ...
    'Rear Rear Swing', 'Rear Rear Lift'};

figure(2); clf; hold on;
for i = 1:numel(act_dof)
    subplot(4,2,i); hold on; title(title_str{i})
    plot(tt_sol, vv_sol_OL(i,:));
    plot(tt_sol, vv_sol_CL(i,:));
    legend('OL Inputs', 'CL Inputs'); 
%     ylim([dp.Vg, dp.Vb])
end

figure(3); clf; hold on;
for i = 1:numel(act_dof)
    subplot(4,2,i); hold on; title(title_str{i})
    plot(tt_sol, xx_sol_OL(act_dof(i), :));
    plot(tt_sol, xx_sol_CL(act_dof(i), :));
    plot(ttd, xxd(act_dof(i), :));
    legend('OL Act Defl', 'CL Act Defl', 'Desired Act Defl');
%     ylim([-SWINGAMP, SWINGAMP])
end

figure(4); clf; hold on;
for i = 1:numel(act_dof)
    subplot(4,2,i); hold on; title(title_str{i})
    plot(tt_sol, xx_sol_OL(nq+act_dof(i), :));
    plot(tt_sol, xx_sol_CL(nq+act_dof(i), :));
    plot(ttd, xxd(nq+act_dof(i), :));
    legend('OL Act Vel', 'CL Act Vel', 'Desired Act Vel');
end

fpopt.loc = 'foot';
fpopt.base_frame = 'World';
pfCL = zeros([numel(tt_sol), size(hamr.HamrRBM.FOOT_POS')]);
vfCL = zeros([numel(tt_sol), size(hamr.HamrRBM.FOOT_POS')]);
nc = numel(hamr.HamrRBM.LEG_NAME);

for i = 1:numel(tt_sol)    
    [pfCL(i,:,:), Ji] = hamr.HamrRBM.getFootPosition(xx_sol_CL(1:nq,i), xx_sol_CL(nq+(1:nv),i), fpopt);
    for j = 1:numel(Ji)
        vfCL(i, :,j) = Ji{j}*xx_sol_CL(nq+(1:nv),i);
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
            plot(tt_sol*1e-3, rad2deg(xx_sol_OL(i,:)), 'LineWidth', 1.5); 
            plot(tt_sol*1e-3, rad2deg(xx_sol_CL(i,:)), 'LineWidth', 1.5); 
        else
            plotyy(tt_sol*1e-3, xx_sol_CL(i,:),tt_sol*1e-3, xx_sol_CL(nq+i,:)); 
        end
    end    
end

tilefigs;

%% Playback 
xtraj_scaled = DTTrajectory(xtraj_sim_CL.getBreaks()*1e-3, xtraj_sim_CL.eval(xtraj_sim_CL.getBreaks()));
xtraj_scaled = xtraj_scaled.setOutputFrame(xtraj_sim_CL.getOutputFrame());
options.slider = true;
v.playback(xtraj_scaled, options);

