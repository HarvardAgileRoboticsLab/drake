clear; clc; close all;

save_dir = '~/Dropbox/CurrentWork/FrictionTrajOpt/MatFiles/TrajOptFiles/';
fname = 'TrajOpt_MovingBody_SimpleSprings10';  

%% Build Simple Robot 

% file
urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'urdf', 'HAMRSimple_scaled.urdf');

% options
options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = true;
options.collision = true;

hamrSimple = HAMRSimpleRBM(urdf,options);
nqS = hamrSimple.getNumPositions; 
nvS = hamrSimple.getNumVelocities; 

%% Load Simple trajectory 

SimpleTraj = load([save_dir, fname]); 
ttSimple = SimpleTraj.xtraj.getBreaks(); 
xxSimple = SimpleTraj.xtraj.eval(ttSimple); 
phiS = zeros(4, numel(ttSimple)); 

for i = 1:numel(ttSimple)
    qS = xxSimple(1:nqS, i);
    qdS = xxSimple(nqS+(1:nvS), i);
    kinsol = hamrSimple.doKinematics(qS, qdS);
    phiS(:,i) = hamrSimple.contactConstraints(kinsol, false);
end

%% Optimize

 [hamr,xtraj,utraj,ctraj,btraj,...
    psitraj,etatraj,jltraj, kltraj, straj, ...
    z,F,info,infeasible_constraint_name] = HAMRVariationalTrajOpt([save_dir, fname], phiS);

save([save_dir, fname, '_Variational'], 'xtraj', 'utraj', 'ctraj', 'btraj', 'psitraj', 'etatraj', ...
    'jltraj', 'kltraj', 'straj'); 

%% Playback

v = hamr.constructVisualizer();
tt = xtraj.getBreaks();
qq = xtraj.eval(tt); qq = qq(1:hamr.getNumPositions(), :); 
qtraj_scaled = PPTrajectory(foh(tt*1e-3, qq));
qtraj_scaled = qtraj_scaled.setOutputFrame(v.getInputFrame());
v.playback(qtraj_scaled, struct('slider', true));

%% Plotting

nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nx = nq+nv;
nu = hamr.getNumInputs();


hh = mean(diff(tt))/2; 
yy = xtraj.eval(tt);
vv = xtraj.eval(tt+hh/2);
uu = utraj.eval(tt);

act_dof = hamr.getActuatedJoints();
nq = hamr.getNumPositions();

figure(1); clf; hold on;
for i = 1:nu
    subplot(4,2,i); hold on; %title(title_str(i))
    yyaxis left; hold on; plot(tt, uu(i,:)); ylabel('Force(N)')
    yyaxis right; hold on; plot(tt, yy(act_dof(i), :)); ylabel('Deflection(mm)')
    legend('Force', 'Deflection')
end

figure(2); clf;
title_str = {'com x', 'com y', 'com z', 'roll', 'pitch', 'yaw'};
for i = 1:6
    subplot(3,2,i); hold on; title(title_str(i))
    yyaxis left; hold on; plot(tt, yy(i,:)); ylabel('Position')
    yyaxis right; hold on; plot(tt+hh/2, vv(i+nq, :)); ylabel('Velocity')
%     legend('Position', 'Velocity')
end

phi = zeros(4, numel(tt));
cc = ctraj.eval(ctraj.getBreaks + hh/2);
for i = 2:numel(tt)
    q = yy(1:nq, i);
    qd = yy(nq+1:nx, i);
    kinsol = hamr.doKinematics(q, qd);
    phi(:,i-1) = hamr.contactConstraints(kinsol, false);
end


legs = {'FLL4'};

figure(3); clf; hold on; 
for i = 1:size(phi, 1)    
    subplot(2,2,i); hold on; %title(legs{i}); 
    yyaxis left; plot(tt, phi(i,:)); ylabel('Distance (mm)'); %ylim([0, 5])
    yyaxis right; plot(tt, cc(i,:)); ylabel('Force (N)')
%     plot(tt, straj.eval(straj.getBreaks()), 'k')
end


figure(4); clf; hold on; 
for i = 1:size(phi, 1)    
    subplot(2,2,i); hold on; %title(legs{i});
    plot(tt, 1e-4*ones(numel(tt),1), 'k'); 
    plot(tt, phi(i,:).*cc(i,:)- straj.eval(straj.getBreaks()));
%     yyaxis right; plot(tt, cc(i,:)); ylabel('Force (N)')
    
end