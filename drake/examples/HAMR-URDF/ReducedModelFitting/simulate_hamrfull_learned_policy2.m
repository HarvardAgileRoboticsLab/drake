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

%% load policy

policy_name = 'ppotraj_new_3';
policy = load(strcat(policy_name, '.mat'));
tp = policy.t;
Vp = double(policy.u');
Vp_hfull = Vp([1, 5, 2, 6, 3, 7, 4, 8],:);

Vtraj = PPTrajectory(zoh(tp, Vp_hfull));
Vtraj = Vtraj.setOutputFrame(hamrWact.getInputFrame);
hamrWact_OL = cascade(Vtraj, hamrWact);


figure(1); clf; hold on;
for i = 1:4
    subplot(4, 2, 2*i-1); hold on;
    title('Swing')
    plot(tp, Vp_hfull(2*i - 1, :));
    subplot(4, 2, 2*i); hold on;
    plot(tp, Vp_hfull(2 * i,:));
    title('Lift')
end


%% Simulate
T = policy.t(end);
tt = 0:options.dt:T;

x0 = zeros(76, 1);
x0(3) = 12.69;
% x0(1:6) = [-6.12311548e-02; 1.36816396e-01; 1.24401108e+01;
%     -2.80173663e-03; -4.19730058e-03; -3.10018578e-03];

[ytraj, xtraj] = simulate(hamrWact_OL, [0, T], x0);

xx = xtraj.eval(tt);
yy = ytraj.eval(tt);

ihip = find(contains(hamr.getStateFrame().getCoordinateNames(), ...
    {'l2', 's2'}) == 1);
xhip = xx(ihip, :);

q = [xx(1:6, :); xhip([1, 3, 5, 7, 2, 4, 6, 8], :)];
qdot = [xx(hamr.getNumPositions() + (1:6), :); ...
    xhip(hamr.getNumInputs() + [1, 3, 5, 7, 2, 4, 6, 8], :) ];

save(strcat(policy_name, '_hfull.mat'), 'q', 'qdot', 'Vp_hfull'); 


%% comparison


qppo = policy.q(:, 1:14)';
qdot_ppo = policy.q(:, 15:end)';


figure(1); clf;
for i = 1:14
    subplot(7, 2, i); hold on;
    title(num2str(i))
    if i == 1
        plot(tt, -q(i, :))
    else
        plot(tt, q(i, :))
    end
        
    plot(tp, qppo(i, :));
    %     legend('
end

figure(2); clf;
for i = 1:14
    subplot(7, 2, i); hold on;
    title(num2str(i))
    if i == 1
        plot(tt, -qdot(i, :))
    else
        plot(tt, qdot(i, :))
    end
    plot(tp, qdot_ppo(i, :));
    %     legend('
end

% figure(3); clf;
% for i = 1:8
%     subplot(4, 2, i); hold on;
%     title(num2str(i))
%     plot(V(i, :))
%     plot(Vp(i, :), '--');
%     %     legend('
% end
% % xx = xtraj_sim_CL.eval(tt);

%% Playback
xtraj_scaled = DTTrajectory(tt*1e-3, xx);
xtraj_scaled = xtraj_scaled.setOutputFrame(xtraj.getOutputFrame());
options.slider = true;
v.playback(xtraj_scaled, options);
%
% %
%
% %% Simulate Open loop
%
%
%
% %% Transmission kin
%
% legs = {'FLL4', 'RLL4', 'FRL4', 'RRL4'};
%
% % Leg positions in leg frame
% leg_pos = [0, 7.58, -11.350;
%     0, 7.58, -11.350;
%     0, -7.58, -11.350;
%     0, -7.58, -11.350];
%
% nc = size(leg_pos,1);
%
% %% Get Actuator Position Velocity
% disp("Calculating actuator pos/vel...")
%
% actuated_dof = hamr.getActuatedJoints();
% qq_act = reshape(qq_full(actuated_dof, :)', N, 2, nc);
% vv_act = reshape(vv_full(actuated_dof, :)', N, 2, nc);
%
% figure(3); clf;
% for i = 1:nc
%     subplot(nc/2, 2, i); hold on;
%     title(legs{i});
%     plot(qq_act(:,1,i), qq_act(:,2,i), '.')
%     xlabel('swing'); ylabel('lift');
%     axis tight;
% end
%
% %% Get hamr simple "hip joint" position velocity
% disp("Calculating hip joint pos/vel...")
%
% ind_simple = find(contains(hamr.getStateFrame().getCoordinateNames(), ...
%     {'l2', 's2'}) == 1);  % index of SFB inputs
%
% qq_simple = reshape(xx(ind_simple(1:8), :)', N, 2, nc);
% vv_simple = reshape(xx(ind_simple(9:16), :)', N, 2, nc);
%
% figure(4); clf;
% for i = 1:nc
%     subplot(nc/2, 2, i); hold on;
%     title(legs{i});
%     plot((180/pi)*qq_simple(:,1,i), (180/pi)*qq_simple(:,2,i), '.')
%     xlabel('swing'); ylabel('lift');
%     axis tight;
% end
%
% %% Get actuator forces
% disp("Calculating actuator forces...")
%
% ffact = reshape(yy(nx+(1:nu), :)', N, 2, nc);
%
% figure(5); clf;
% for i = 1:nc
%     subplot(nc/2, 2, i); hold on;
%     title(legs{i});
%     plot(ffact(:,1,i), ffact(:,2,i), '.')
%     xlabel('swing'); ylabel('lift');
%     axis tight;
% end
%
% %% Get equivalent hip torques
% disp("Calculating hip torques...")
%
% tthip = zeros(2*nc, N);
% J_sfb_a = zeros(nq, numel(actuated_dof));
% J_sfb_a(actuated_dof, :) = eye(numel(actuated_dof));
%
% hamr_state_frame = hamr.getStateFrame();
% unactuated_dof = find(contains(hamr_state_frame.getFrameByName(...
%     'HamrPosition').getCoordinateNames(), 'act') == 0);  % index of unactuated dofs
%
% for j = 1:N
%     [K, dK] = hamr.positionConstraints(xx(1:nq, j));
%     dKUnique = dK(hamr.valid_loops, :);
%     Jc = -dKUnique(:,unactuated_dof)\dKUnique(:, actuated_dof);
%     J_sfb_a(unactuated_dof, :) = Jc;
%     tthip(:, j) = J_sfb_a(ind_simple(1:8), :)' * yy(nx+(1:nu), j);
% end
%
% tthip = reshape(tthip', N, 2, nc);
%
% figure(6); clf;
% for i = 1:nc
%     subplot(nc/2, 2, i); hold on;
%     title(legs{i});
%     plot(tthip(:,1,i), tthip(:,2,i), '.')
%     xlabel('swing'); ylabel('lift');
%     axis tight;
% end
%
% %% Get Leg Position and Velocity
% disp("Calculating leg pos/vel...")
%
% fkopt.base_or_frame_id = hamr.findLinkId('Chassis');
%
% qq_leg = zeros([numel(tt), size(leg_pos')]);
% vv_leg = zeros([numel(tt), size(leg_pos')]);
%
% % Compute leg positions
% for j = 1:N
%
%     kinsol = hamr.doKinematics(xx(1:nq, j), xx(nq+(1:nv), j));
%     for i = 1:size(leg_pos,1)
%         [qq_leg(j,:,i), J] = hamr.forwardKin(kinsol, hamr.findLinkId(legs{i}), ...
%             leg_pos(i,:)', fkopt);
%         vv_leg(j, :, i) = J* xx(nq+(1:nv), j);
%     end
%
% end
%
% figure(7); clf;
% nc = size(leg_pos,1);
% for i = 1:nc
%
%     subplot(nc/2, 2, i); hold on;
%     title(legs{i}); view(3);
%     plot3(qq_leg(:,1,i), qq_leg(:,2,i), qq_leg(:,3,i), '.')
%     xlabel('X'); ylabel('Y'); zlabel('Z')
%     axis tight;
%
% end
%
% %% Save
% disp("saving data ...")
%
% data = struct();
%
% data.time = tt;
%
% % Inputs
% data.inputs.voltage = reshape(VV', N, 2, nc);
% data.inputs.actuator_force = ffact;
% data.inputs.hip_torque = tthip;
%
% % state
% if options.floating
%     data.state.qfloat = qq_full(1:6, :);
%     data.state.vfloat = vv_full(1:6, :);
%     data.state.qfull = reshape(qq_full(7:end, :)', N, 8, nc);
%     data.state.vfull = reshape(vv_full(7:end, :)', N, 8, nc);
% else
%     data.state.qfull = reshape(qq_full', N, 8, nc);
%     data.state.vfull = reshape(vv_full', N, 8, nc);
% end
% data.state.qact = qq_act;
% data.state.vact = vv_act;
% data.state.qhip = qq_simple;
% data.state.vhip = vv_simple;
%
% % outputs
% data.outputs.qleg = qq_leg;
% data.outputs.vleg = vv_leg;
%
% save('hamr_data_pronk2', 'data')
%
% %% Playback
% xtraj_scaled = DTTrajectory(xtraj.getBreaks()*1e-3, xtraj.eval(xtraj.getBreaks()));
% xtraj_scaled = xtraj_scaled.setOutputFrame(xtraj.getOutputFrame());
% options.slider = true;
% v.playback(xtraj_scaled, options);
%
%
%
