clear; clc; close all; 
save_dir = '~/Dropbox/CurrentWork/FrictionTrajOpt/MatFiles/StepUp/';
fname = 'StepUp_5_mm'; 

trajVar = load([save_dir, fname, '.mat']); 

ttVar = trajVar.xtraj.getBreaks(); 
hhVar = mean(diff(ttVar)); 
uuVar = trajVar.utraj.eval(ttVar); 
xxVar = trajVar.xtraj.eval(ttVar + hhVar/2); 

%% Box
w  = 100; 
l = 50;
h = 5; 
boxes = [0, 0, l, w, h]; 

%% Build Full Model

urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2.urdf');

options.terrain = RigidBodyStepTerrain(boxes);
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = true; %false;
options.collision = true; %false;
options.dt = 1; %0.1; 

% Build robot + visualizer
hamr = HamrTSRBM(urdf, options);
hamr = compile(hamr);
nq = hamr.getNumPositions(); 
nu = hamr.getNumInputs(); 
v = hamr.constructVisualizer(); 
act_dof = hamr.getActuatedJoints(); 

%% Plotting 

figure(1); clf; hold on;
for i = 1:nu
    subplot(4,2,i); hold on; %title(title_str(i))'
    plot(ttVar+hhVar/2, uuVar(i,:)); ylabel('Force(N)')
%     plot(ttID+hhID/2, uuID(i, :)); 
end
% 
figure(2); clf; hold on;
for i = 1:nu
    subplot(4,2,i); hold on; %title(title_str(i))'
    plot(ttVar, xxVar(act_dof(i),:)); ylabel('Deflection')
%     plot(ttID, xxID(act_dof(i)-6, :)); 
end

%% Simulate
xtraj_playback = PPTrajectory(foh(ttVar/1e3, xxVar(1:nq, :))); 
xtraj_playback = xtraj_playback.setOutputFrame(v.getInputFrame); 
v.playback(xtraj_playback, struct('slider', true))
