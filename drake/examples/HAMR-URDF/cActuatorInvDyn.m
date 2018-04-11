clear; clc; close all;

%% Build Robot

urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2_TYM.urdf');

% options
options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = true;
options.collision = true;
options.z_inactive_guess_tol = 0.1;
options.dt = 1;

% Build robot + visualizer
hamr = HamrTSRBM(urdf, options);
hamr = compile(hamr);

%% Trajectory
save_dir = '~/Dropbox/CurrentWork/FrictionTrajOpt/MatFiles/SimWarmStart/';
fname = 'TROT_0.1N_30Hz_TYM_TEF_VariationalSmooth_converted';
append = '';
trajTrans = load([save_dir, fname, append, '.mat']); %, '_VariationalMU.mat']);

% Time
tt = trajTrans.xtraj.getBreaks();
hh = mean(diff(tt));

% State
xtraj = trajTrans.xtraj;
xx = xtraj.eval(tt);
act_dof = hamr.getActuatedJoints();
xx_act = xx(act_dof, :);
xx_act_m = (xx_act(:, 1:end-1) + xx_act(:, 2:end))/2; 

% input
utraj = trajTrans.utraj;
uu = utraj.eval(tt(1:end-1));

% Rest
ctraj = trajTrans.ctraj;
btraj = trajTrans.btraj;
psitraj = trajTrans.psitraj;
etatraj = trajTrans.etatraj;
jltraj = trajTrans.jltraj;
kltraj = trajTrans.kltraj;
straj = trajTrans.straj;

%% Build Actuators

dp.Vb = 225;
dp.Vg = 0;

nact = 8;
act_names = {'FLsact', 'FLlact', 'RLsact', 'RLlact', ...
    'FRsact', 'FRlact', 'RRsact', 'RRlact'};
hr_actuators = HamrActuators(nact, act_names,  [-1; -1; 1; 1; -1; -1; 1; 1], dp);

for i = 2:2:nact
    hr_actuators.dummy_bender(i) = hr_actuators.dummy_bender(i).setCFThickness(0.1); 
end

%%
vv = zeros(size(uu));
for i = 1:(numel(tt)-1)
    [vv(:,i), dvv] = ipzt_fun(hamr, hr_actuators, tt(i), xx_act_m(:,i), uu(1:(nu-nl),i));
    vv(vv(:,i) < 0, i) = 0; 
    vv(vv(:,i) > 225, i) = 225; 
end

figure(1); clf; hold on;
for i = 1:nact
    subplot(nact/2, 2, i)
    plot(tt(1:end-1), vv(i,:)); title(act_names{i})
    ylim([0, 225])
end

vtraj = PPTrajectory(zoh(tt(1:end-1)+hh/2, vv));
save([save_dir, fname, append, '_PlusAct.mat'], 'xtraj', 'utraj', 'vtraj', ...
    'ctraj', 'btraj', 'psitraj', 'etatraj', 'jltraj', 'kltraj', 'straj');




