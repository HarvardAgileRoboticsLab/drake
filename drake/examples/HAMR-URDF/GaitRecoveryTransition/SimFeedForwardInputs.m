clear; clc; close all;
warning('off', 'MATLAB:nargchk:deprecated');
addpath('../', '../urdf/')
DL = linspace(-25, 75, 5); 

%% Load Data

data_dir = 'data_16-May-2018';
freq = 10;

ffinputs = load([data_dir, '/FFInputs_', num2str(freq), 'Hz.mat']);
ffinputs = ffinputs.ffinputs;

N = randperm(size(ffinputs, 1), 1) %size(ff_inputs, 1);
M = randperm(size(ffinputs, 2), 1) %size(ff_inputs, 2);
TYPE = 2; 

ff_inputs_nm = ffinputs{N, M};


for i = 1:size(ffinputs, 1)
    for j = 1:size(ffinputs, 2)
%     ffinputs{i,j}.params.DL = DL(j);
%     DLmat(i,j) =  ffinputs{i,j}.params.DL;
    end
end

% save([data_dir, '/FFInputs_', num2str(freq), 'Hz.mat'], 'ffinputs')
%% Load Rigid Body Model

urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2.urdf');

% options
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.z_inactive_guess_tol = 0.1;
options.use_bullet = false;

% options to change
dt = 0.4;
options.dt = dt;
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
nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nu = hamr.getNumInputs();
qa = hamr.getActuatedJoints();
nqa = numel(qa); 

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

%% Extract Inputs

tt_sol = ff_inputs_nm.ttopt;
vv_sol = ff_inputs_nm.vvopt; 
xx_sol = ff_inputs_nm.xxopt;
% NPTS = numel(tt_sol); 
xx_des = GenerateHueristicActTraj(ff_inputs_nm.params.GAIT, freq, ...
    ff_inputs_nm.params.DC, ff_inputs_nm.params.DL, ff_inputs_nm.params.NSAMP, TYPE);

xx_des(:, 2:2:end) = ff_inputs_nm.params.SWINGAMP*xx_des(:, 2:2:end);
xx_des(:, 3:2:end) = ff_inputs_nm.params.LIFTAMP*xx_des(:, 3:2:end);

%% Plot and Simulate 
NCYC = 5; 
tt_plot = linspace(0, NCYC*tt_sol(end), NCYC*numel(tt_sol));
vv_plot = repmat(vv_sol, 1, NCYC);

figure(1); clf;
for l = 1:nqa
    subplot(nqa/2, 2, l); hold on;
    plot(tt_sol, xx_sol(qa(l),:))
    plot(1e3*xx_des(:,1), xx_des(:,l+1), 'k')
%     plot(linspace(0, in_params.T, NSAMP), trajd(:,l+1)', 'k--')
end
drawnow;

figure(2); clf;
NCYC = 2;
for l = 1:nqa
    subplot(nqa/2, 2, l); hold on;
    %     plotyy(ttopt, uuopt(j,:), ttopt, vvopt(j,:))
    plot(tt_plot, vv_plot(l,:), '-')
    ylim([hr_actuators.dummy_bender(1).dp.Vg, hr_actuators.dummy_bender(1).dp.Vb])
end
drawnow;

%% Simulate

tramp = 2/(freq/1e3);
ramp = tt_plot/tramp; ramp(tt_plot >= tramp) = 1;

vv_sim =  ramp.*(vv_plot - hr_actuators.dummy_bender(1).dp.Vb/2) + ...
    hr_actuators.dummy_bender(1).dp.Vb/2;

vv_sim(1:2:end, :) =  -(vv_sim(1:2:end, :) - hr_actuators.dummy_bender(1).dp.Vb/2) + ...
    hr_actuators.dummy_bender(1).dp.Vb/2;

vtraj = PPTrajectory(foh(tt_plot, vv_sim));
vtraj = setOutputFrame(vtraj, hamrWact.getInputFrame());

hamr_OL = cascade(vtraj, hamrWact);
xtraj_sim = simulate(hamr_OL, [0 tt_plot(end)], x0);

%% Playback
xtraj_scaled = DTTrajectory(xtraj_sim.getBreaks()*1e-3, xtraj_sim.eval(xtraj_sim.getBreaks()));
xtraj_scaled = xtraj_scaled.setOutputFrame(xtraj_sim.getOutputFrame());
options.slider = true;
v.playback(xtraj_scaled, options);




