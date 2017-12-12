clear; clc; close all; 
save_dir = '~/Dropbox/CurrentWork/FrictionTrajOpt/MatFiles/TrajOptFiles/';
fname = 'TrajOpt_MovingBody_SimpleSprings10'; 

SimpleTraj = load([save_dir, fname]); 
IDTraj = load([save_dir, fname, '_fullRobot']); 
VarTraj = load([save_dir, fname, '_Variational']); 

ttSimple = SimpleTraj.xtraj.getBreaks(); 
hhSimple = mean(diff(ttSimple)); 
uuSimple = SimpleTraj.utraj.eval(ttSimple); 

ttID = IDTraj.xtraj.FL_scaled.getBreaks(); 
hhID = mean(diff(ttID)); 

xtrajID = [IDTraj.xtraj.FL_scaled; 
    IDTraj.xtraj.RL_scaled; 
    IDTraj.xtraj.FR_scaled; 
    IDTraj.xtraj.RR_scaled];
xxID = xtrajID.eval(ttID); 

utrajID = [IDTraj.utraj.FL_scaled(1:2);
    IDTraj.utraj.RL_scaled(1:2);
    IDTraj.utraj.FR_scaled(1:2);
    IDTraj.utraj.RR_scaled(1:2)];
uuID = utrajID.eval(ttID+hhID/2); 

ttVar = VarTraj.xtraj.getBreaks(); 
hhVar = mean(diff(ttVar)); 
uuVar = VarTraj.utraj.eval(ttVar); 
xxVar = VarTraj.xtraj.eval(ttVar); 

%% Build Full Model

urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2.urdf');

options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = true; %false;
options.collision = true; %false;
options.dt = 0.5; %0.1; 

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
    plot(ttSimple+hhSimple/2, uuSimple(i,:));
    plot(ttVar+hhVar/2, uuVar(i,:)); ylabel('Force(N)')
    plot(ttID+hhID/2, uuID(i, :)); 
end
% 
figure(2); clf; hold on;
for i = 1:nu
    subplot(4,2,i); hold on; %title(title_str(i))'
    plot(ttVar, xxVar(act_dof(i),:)); ylabel('Deflection')
    plot(ttID, xxID(act_dof(i)-6, :)); 
end