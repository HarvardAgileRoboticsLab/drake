clear; clc; close all; 
save_dir = '~/Dropbox/CurrentWork/FrictionTrajOpt/MatFiles/AgileBehaviors/';
fname = 'SINGLE_JUMP_0.25N_2Hz_TYM'; 

OrigTraj = load([save_dir, fname, '.mat']); 
VarTraj = load([save_dir, fname, 'b.mat']); 

ttSimple = OrigTraj.tt;
hhSimple = mean(diff(ttSimple)); 
xxSimple = OrigTraj.xx; 
uuSimple = OrigTraj.uu; 

% last cycle
tokens = strsplit(fname, '_'); 
freq = str2double(tokens{3}(1:end-2)); 
nc = ttSimple(end)*freq*1e-3; 
nCF0 = find(ttSimple >= (nc-1)/freq/1e-3, 1, 'first');
nCF1 = numel(ttSimple); 

ttSimple = ttSimple(nCF0:nCF1) - ttSimple(nCF0);
uuSimple = uuSimple(:, nCF0:nCF1);
xxSimple = xxSimple(:, nCF0:nCF1);

% zero out starting x and y pos
xxSimple0 = 0*xxSimple(:,1); xxSimple0(1:2) = xxSimple(1:2, 1); 
xxSimple = bsxfun(@minus, xxSimple, xxSimple0); 

ttVar = VarTraj.xtraj.getBreaks(); 
hhVar = mean(diff(ttVar)); 
uuVar = VarTraj.utraj.eval(ttVar); 
xxVar = VarTraj.xtraj.eval(ttVar + hhVar/2); 

%% Build Full Model

urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2.urdf');

options.terrain = RigidBodyFlatTerrain();
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
nv = hamr.getNumVelocities(); 
nu = hamr.getNumInputs(); 
v = hamr.constructVisualizer(); 
act_dof = hamr.getActuatedJoints(); 

%% Plotting 

figure(1); clf; hold on;
for i = 1:nu
    subplot(4,2,i); hold on; %title(title_str(i))'
    plot(ttSimple, uuSimple(i,:));
    plot(ttVar+hhVar/2, uuVar(i,:)); ylabel('Force(N)')
%     plot(ttID+hhID/2, uuID(i, :)); 
end
% 
figure(2); clf; hold on;
for i = 1:nu
    subplot(4,2,i); hold on; %title(title_str(i))'
    plot(ttSimple, xxSimple(act_dof(i), :))
    plot(ttVar, xxVar(act_dof(i),:)); ylabel('Deflection')
%     plot(ttID, xxID(act_dof(i)-6, :)); 
end


phi = zeros(4, numel(ttVar));
cc = VarTraj.ctraj.eval(ttVar + hhVar/2);
for i = 2:numel(ttVar)
    q = xxVar(1:nq, i);
    qd = xxVar(nq+(1:nv), i);
    kinsol = hamr.doKinematics(q, qd);
    phi(:,i-1) = hamr.contactConstraints(kinsol, false);
end

figure(3); clf; hold on; 
for i = 1:size(phi, 1)    
    subplot(2,2,i); hold on; %title(legs{i}); 
    yyaxis left; plot(ttVar, phi(i,:)); ylabel('Distance (mm)'); %ylim([0, 5])
    yyaxis right; plot(ttVar, cc(i,:)); ylabel('Force (N)')
%     plot(tt, straj.eval(straj.getBreaks()), 'k')
end

xtraj_playback = PPTrajectory(foh(ttVar/1e3, xxVar(1:nq, :)));
xtraj_playback = xtraj_playback.setOutputFrame(v.getInputFrame()); 
v.playback(xtraj_playback, struct('slider', true))