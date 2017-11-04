clear; clc; close all;

%%

urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2.urdf');

% options
options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = false;
options.collision = false;
options.z_inactive_guess_tol = 0.1;
options.dt = 1; 

% Build robot + visualizer
hamr = HamrTSRBM(urdf, options);
hamr = compile(hamr);

%% Trajectory
fname = 'TrajOpt_26-Oct-2017_6_fullRobot'; 
traj_full = load(fname); 

% Build transmission trajectory
xtraj = traj_full.xtraj();
tt = xtraj.FL_scaled.getBreaks();
hh = mean(diff(tt));

xxFL = xtraj.FL_scaled.eval(tt);
xxRL = xtraj.RL_scaled.eval(tt);
xxFR = xtraj.FR_scaled.eval(tt);
xxRR = xtraj.RR_scaled.eval(tt);

% Build input trajectory (included contact and const. forces)
utraj = traj_full.utraj;
uFL = utraj.FL_scaled.eval(tt+hh/2);
uRL = utraj.RL_scaled.eval(tt+hh/2);
uFR = utraj.FR_scaled.eval(tt+hh/2);
uRR = utraj.RR_scaled.eval(tt+hh/2);
nut = size(utraj.FL_scaled.eval(tt), 1);

%% Just inputs and actuated dof

act_dof = [1,5]; 
xx = [xxFL(act_dof, :); xxRL(act_dof, :); xxFR(act_dof, :); xxRR(act_dof, :)]; 
uu = [uFL(1:2, :); uRL(1:2,:); uFR(1:2,:); uRR(1:2, :)];

%% Build Actuators 

dp.Vb = 300;
dp.Vg = 0;

nact = 8; 
hr_actuators = HamrActuators(nact, {'FLsact', 'FLlact', 'RLsact', 'RLlact', ...
    'FRsact', 'FRlact', 'RRsact', 'RRlact'},  [1; 1; -1; -1; 1; 1; -1; -1], dp);

%% 
vv = zeros(size(uu)); 
for i = 1:numel(tt)
    [vv(:,i), dvv] = ipzt_fun(hamr, hr_actuators, tt(i), xx(:,i), uu(:,i));
end

vtraj = PPTrajectory(zoh(tt+hh/2, vv)); 
save([fname, 'PlusAct'], 'xtraj', 'utraj', 'vtraj'); 
