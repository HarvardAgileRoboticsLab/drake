clear; clc; close all; 

load_dir = '/media/nddoshi/NeelWoodlab1/FrictionTrajOpt/JumpData/'; 
file = 'SINGLE_JUMP_0.25N_10Hz_TYM';
traj = load([load_dir, file, '.mat']);

%% Build Robot
urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2_TYM.urdf');

options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = true; %false;
options.collision = true; %false;
options.dt = 1; %0.1;

% Build robot + visualizer
hamr = HamrTSRBM(urdf, options);
nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nu = hamr.getNumInputs();
nc = hamr.getNumContactPairs();
v = hamr.constructVisualizer();

%% Make trajectory
NO = numel(traj.tt); 
NCYC = 1; 
T = NCYC*traj.tt(end);
tt = 1e-3*linspace(traj.tt(1), T, NCYC*numel(traj.tt));
xx = repmat(traj.xx, 1, NCYC); 

for i = 2:NCYC
    xx(1,(i-1)*NO+1:i*NO) = xx(1, (i-1)*NO+1:i*NO) + xx(1,(i-1)*NO);
end

playback_traj = PPTrajectory(foh(tt, xx(1:nq,:))); 
playback_traj = playback_traj.setOutputFrame(v.getInputFrame); 

%% Make visualizer
v.playback(playback_traj, struct('slider', true)); 



