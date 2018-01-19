clear; clc; close all;

%% Build Robot

urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2.urdf');

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
fname = 'TROT_0.2N_10Hz_MU10';
append = '_VariationalSmooth2';
trajTrans = load([save_dir, fname, append, '.mat']); %, '_VariationalMU.mat']);

% Time
tt = trajTrans.xtraj.getBreaks();
hh = mean(diff(tt));

% State
xtraj = trajTrans.xtraj;
xx = xtraj.eval(tt);
act_dof = hamr.getActuatedJoints();
xx_act = xx(act_dof, :);

% input
utraj = trajTrans.utraj;
uu = utraj.eval(tt);

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

%%
vv = zeros(size(uu));
for i = 1:numel(tt)
    [vv(:,i), dvv] = ipzt_fun(hamr, hr_actuators, tt(i), xx_act(:,i), uu(:,i));
end

figure(1); clf; hold on;
for i = 1:nact
    subplot(nact/2, 2, i)
    plot(tt, vv(i,:)); title(act_names{i})
    ylim([0, 225])
end

vtraj = PPTrajectory(zoh(tt, vv));
save([save_dir, fname, append, '_PlusAct.mat'], 'xtraj', 'utraj', 'vtraj', ...
    'ctraj', 'btraj', 'psitraj', 'etatraj', 'jltraj', 'kltraj', 'straj');




