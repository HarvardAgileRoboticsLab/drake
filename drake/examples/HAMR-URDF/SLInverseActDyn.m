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
fname = 'TROT_0.2N_10Hz';
trajTrans = load([save_dir, fname, '_Variational.mat']); %, '_VariationalMU.mat']);

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
uu = utraj.eval(tt + hh/2); 

% Rest
ctraj = trajTrans.ctraj;
btraj = trajTrans.btraj; 
psitraj = trajTrans.psitraj; 
etatraj = trajTrans.etatraj; 
jltraj = trajTrans.jltraj; 
kltraj = trajTrans.kltraj; 
straj = trajTrans.straj; 

%% Build Actuators 

dp.Vb = 200;
dp.Vg = 0;

nact = 8; 
hr_actuators = HamrActuators(nact, {'FLsact', 'FLlact', 'RLsact', 'RLlact', ...
    'FRsact', 'FRlact', 'RRsact', 'RRlact'},  [1; 1; -1; -1; 1; 1; -1; -1], dp);

%% 
vv = zeros(size(uu)); 
for i = 1:numel(tt)
    [vv(:,i), dvv] = ipzt_fun(hamr, hr_actuators, tt(i), xx_act(:,i), uu(:,i));
end

vtraj = PPTrajectory(zoh(tt, vv)); 
save([save_dir, fname, '_VariationalPlusAct.mat'], 'xtraj', 'utraj', 'vtraj', ...
    'ctraj', 'btraj', 'psitraj', 'etatraj', 'jltraj', 'kltraj', 'straj'); 
