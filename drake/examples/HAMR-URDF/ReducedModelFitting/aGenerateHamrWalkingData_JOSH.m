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
gait = 'TROT';
ISFLOAT = true; % floating (gnd contact) or in air (not floating)
options.floating = ISFLOAT;
options.collision = ISFLOAT;
x0 = zeros(76, 1); x0(3) = 12.69;
options.terrain = RigidBodyFlatTerrain();


% Build robot + visualizer
hamr = HamrTSRBM(HamrRBM(urdf, options), options.dt, options);
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
fd = 0.060;         % drive frequency (Hz)
T = 10/fd;
tt = 0:options.dt:T;
Vmax = 50;
N = numel(tt);

switch gait
    case 'TROT'
        VV = Vmax*[sin(2*pi*fd*tt + 3*pi/2);            % FLswing
            sin(2*pi*fd*tt + pi);                       % FLlift
            sin(2*pi*fd*tt + pi/2);              % RLSwing
            sin(2*pi*fd*tt + pi);                       % RLLift
            sin(2*pi*fd*tt + 3*pi/2);                % FRswing
            sin(2*pi*fd*tt + pi);                       % FRlift
            sin(2*pi*fd*tt + pi/2);                % RRSwing
            sin(2*pi*fd*tt + pi)];               % RRLift
    case 'PRONK'
        VV = Vmax*[sin(2*pi*fd*tt + 3*pi/2);            % FLswing
            sin(2*pi*fd*tt + pi);                       % FLlift
            sin(2*pi*fd*tt - pi/2);              % RLSwing
            sin(2*pi*fd*tt);                       % RLLift
            sin(2*pi*fd*tt + pi/2);              % FRswing
            sin(2*pi*fd*tt);                    % FRlift
            sin(2*pi*fd*tt + pi/2);              % RRSwing
            sin(2*pi*fd*tt + pi)];                      % RRLift
    otherwise
        VV = zeros(8, numel(tt));
end

% ramp
tramp = 2/fd;
ramp = tt/tramp; ramp(tt >= tramp) = 1;

VV = bsxfun(@times, ramp, VV) + 0.5*(dp.Vb - dp.Vg);
Vtraj = PPTrajectory(foh(tt, VV));
Vtraj = setOutputFrame(Vtraj, hamrWact.getInputFrame());

figure(1); clf;
plot(tt, VV(1,:), tt, VV(2,:), '--');
plot(tt, VV(3,:), tt, VV(4,:), '--');

legend('Swing Drive', 'Lift Drive')

%% Simulate Open loop

relative_speed = zeros(10,1);

for i = 1:10
    hamrWact_OL = cascade(Vtraj, hamrWact);

    disp('Simulating...');
    tic; [ytraj, xtraj] = simulate(hamrWact_OL, [0 T], x0); tlcp = toc;
    fprintf('Sim time was %f, \n', tlcp)
    relative_speed(i) = (tlcp/T)*1000;

    % extract
    xx = xtraj.eval(tt);
    yy = ytraj.eval(tt);

    % full transmission state
    qq_full = xx(1:nq, :);
    vv_full = xx(nq+(1:nv), :);
end
%%

figure(100); hold on; 
plot(relative_speed, '*')
plot(1:10, mean(relative_speed)*ones(10,1), 'k')


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
    data.state.qfloat = qq_full(1:6, :);
    data.state.vfloat = vv_full(1:6, :);
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

save('hamr_data_pronk2', 'data')

%% Playback
xtraj_scaled = DTTrajectory(xtraj.getBreaks()*1e-3, xtraj.eval(xtraj.getBreaks()));
xtraj_scaled = xtraj_scaled.setOutputFrame(xtraj.getOutputFrame());
options.slider = true;
v.playback(xtraj_scaled, options);



