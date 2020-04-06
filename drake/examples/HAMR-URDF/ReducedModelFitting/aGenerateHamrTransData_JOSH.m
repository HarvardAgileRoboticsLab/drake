clear; clc; close all;
addpath('../')
% save_dir = '~/Dropbox/CurrentWork/FrictionTrajOpt/MatFiles/SysIDFiles/';

global u_traj kl_traj
u_traj = [];
kl_traj = [];

% Sample Rate and Time
dt = 0.25;
Fs = 1/dt;
T = 8000;
tt = 0:dt:T;
N = numel(tt);

% White Noise Freq Characteristics
f1 = 40e-3;
band = [0, 2*f1*dt];

%% Build robot

% options
name = 'HAMR_scaledV2';
urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'urdf', [name, '.urdf']);
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = false;
options.z_inactive_guess_tol = 0.1;
options.dt = dt;

% build hamr
hamr = HamrTSRBM(urdf, options);
hamr = compile(hamr);
v = hamr.constructVisualizer();

% state/input dimenisons
nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nx = nq+nv;
nu = hamr.getNumInputs();

% Relax joint limits
hamr = hamr.setJointLimits(-Inf(nq, 1), Inf(nq, 1));
hamr = compile(hamr);

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


%% Build hamr with actuators

hamrWact = mimoFeedback(hr_actuators, hamr, connection1, connection2, ...
    input_select, output_select);

%% Build and Plot Input

tramp = 500;
Vscale = (dp.Vb - dp.Vg)/3.5;
VV = Vscale*idinput([N, hamr.getNumInputs()], 'rgs', band)';
VV(tt < 500) = (tt(tt<500)/tramp).*VV(tt<500);
VV = VV + (dp.Vb - dp.Vg)/2 + dp.Vg;
VV(VV > dp.Vb) = dp.Vb;
VV(VV < dp.Vg) = dp.Vg;


figure(1); clf;
subplot(2,1,1);
yyaxis left; hold on; plot(tt, VV(1,:))
subplot(2,1,2);
yyaxis left;  hold on; plot(tt, VV(2,:))

freq = [-(Fs/2):(1/N/dt):(Fs/2)];

f0 = 0;
ind0 = find(freq > f0, 1, 'first');
ind1 = find(freq > 2*f1, 1, 'first');
Vfft = abs(fftshift((fft(bsxfun(@minus, VV, (dp.Vb - dp.Vg)/2 + dp.Vg), [], 2)/size(VV,2))));

figure(2); clf;
subplot(2,1,1);
yyaxis left; hold on; plot(freq(ind0:ind1)*1e3, 2*Vfft(1,ind0:ind1))
subplot(2,1,2); hold on;
yyaxis left; hold on; plot(freq(ind0:ind1)*1e3, 2*Vfft(2,ind0:ind1))

%% IC and  Simulate

Vtraj = PPTrajectory(zoh(tt, VV));
Vtraj = Vtraj.setOutputFrame(hamrWact.getInputFrame);
hamrWact_OL = cascade(Vtraj, hamrWact);

q0 = zeros(nq,1);
x0 = [q0; 0*q0];

disp('Simulating...');
tic; [ytraj, xtraj] = simulate(hamrWact_OL, [0 T], x0); tlcp = toc;
fprintf('Sim time was %f, \n', tlcp)

% extract
xx = xtraj.eval(tt);
yy = ytraj.eval(tt);

% full transmission state
qq_full = xx(1:nq, :);
vv_full = xx(nq+(1:nv), :);

%% Transmission kin

legs = {'FLL4', 'RLL4', 'FRL4', 'RRL4'};

% Leg positions in leg frame
leg_pos = [0, 7.58, -11.350;
    0, 7.58, -11.350;
    0, -7.58, -11.350;
    0, -7.58, -11.350];

nc = size(leg_pos,1);

%% Get Actuator Position Velocity
disp("Calculating actuator pos/vel...")

actuated_dof = hamr.getActuatedJoints();
qq_act = reshape(qq_full(actuated_dof, :)', N, 2, nc);
vv_act = reshape(vv_full(actuated_dof, :)', N, 2, nc);

figure(3); clf;
for i = 1:nc
    subplot(nc/2, 2, i); hold on;
    title(legs{i});
    plot(qq_act(:,1,i), qq_act(:,2,i), '.')
    xlabel('swing'); ylabel('lift');
    axis tight;
end

%% Get hamr simple "hip joint" position velocity
disp("Calculating hip joint pos/vel...")

ind_simple = find(contains(hamr.getStateFrame().getCoordinateNames(), ...
    {'l2', 's2'}) == 1);  % index of SFB inputs

qq_simple = reshape(xx(ind_simple(1:8), :)', N, 2, nc);
vv_simple = reshape(xx(ind_simple(9:16), :)', N, 2, nc);

figure(4); clf;
for i = 1:nc
    subplot(nc/2, 2, i); hold on;
    title(legs{i});
    plot((180/pi)*qq_simple(:,1,i), (180/pi)*qq_simple(:,2,i), '.')
    xlabel('swing'); ylabel('lift');
    axis tight;
end

%% Get actuator forces
disp("Calculating actuator forces...")

ffact = reshape(yy(nx+(1:nu), :)', N, 2, nc);

figure(5); clf;
for i = 1:nc
    subplot(nc/2, 2, i); hold on;
    title(legs{i});
    plot(ffact(:,1,i), ffact(:,2,i), '.')
    xlabel('swing'); ylabel('lift');
    axis tight;
end

%% Get equivalent hip torques
disp("Calculating hip torques...")

tthip = zeros(2*nc, N);
J_sfb_a = zeros(nq, numel(actuated_dof));
J_sfb_a(actuated_dof, :) = eye(numel(actuated_dof));

hamr_state_frame = hamr.getStateFrame();
unactuated_dof = find(contains(hamr_state_frame.getFrameByName(...
    'HamrPosition').getCoordinateNames(), 'act') == 0);  % index of unactuated dofs

for j = 1:N
    [K, dK] = hamr.positionConstraints(xx(1:nq, j));
    dKUnique = dK(hamr.valid_loops, :);
    Jc = -dKUnique(:,unactuated_dof)\dKUnique(:, actuated_dof);
    J_sfb_a(unactuated_dof, :) = Jc;
    tthip(:, j) = J_sfb_a(ind_simple(1:8), :)' * yy(nx+(1:nu), j);
end

tthip = reshape(tthip', N, 2, nc);

figure(6); clf;
for i = 1:nc
    subplot(nc/2, 2, i); hold on;
    title(legs{i});
    plot(tthip(:,1,i), tthip(:,2,i), '.')
    xlabel('swing'); ylabel('lift');
    axis tight;
end

%% Get Leg Position and Velocity
disp("Calculating leg pos/vel...")

fkopt.base_or_frame_id = hamr.findLinkId('Chassis');

qq_leg = zeros([numel(tt), size(leg_pos')]);
vv_leg = zeros([numel(tt), size(leg_pos')]);

% Compute leg positions
for j = 1:N
    
    kinsol = hamr.doKinematics(xx(1:nq, j), xx(nq+(1:nv), j));
    for i = 1:size(leg_pos,1)
        [qq_leg(j,:,i), J] = hamr.forwardKin(kinsol, hamr.findLinkId(legs{i}), ...
            leg_pos(i,:)', fkopt);
        vv_leg(j, :, i) = J* xx(nq+(1:nv), j);
    end
    
end

figure(7); clf;
nc = size(leg_pos,1);
for i = 1:nc
    
    subplot(nc/2, 2, i); hold on;
    title(legs{i}); view(3);
    plot3(qq_leg(:,1,i), qq_leg(:,2,i), qq_leg(:,3,i), '.')
    xlabel('X'); ylabel('Y'); zlabel('Z')
    axis tight;
    
end

%% Save
disp("saving data ...")

data = struct();

data.time = tt;

% Inputs
data.inputs.voltage = reshape(VV', N, 2, nc);
data.inputs.actuator_force = ffact;
data.inputs.hip_torque = tthip;

% state
if options.floating
    data.state.qfloat = qq_full(1:6, :)
    data.state.vfloat = vv_full(1:6, :)
    data.state.qfull = reshape(qq_full(7:end, :)', N, 8, nc);
    data.state.vfull = reshape(vv_full(7:end, :)', N, 8, nc);
else
    data.state.qfull = reshape(qq_full', N, 8, nc);
    data.state.vfull = reshape(vv_full', N, 8, nc);
end
data.state.qact = qq_act;
data.state.vact = vv_act;
data.state.qhip = qq_simple;
data.state.vhip = vv_simple;

% outputs
data.outputs.qleg = qq_leg;
data.outputs.vleg = vv_leg;

save('hamr_data_0grav_fixed', 'data')

%% Playback
xtraj_scaled = DTTrajectory(xtraj.getBreaks()*1e-3, xtraj.eval(xtraj.getBreaks()));
xtraj_scaled = xtraj_scaled.setOutputFrame(xtraj.getOutputFrame());
options.slider = true;
v.playback(xtraj_scaled, options);



