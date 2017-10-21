clear; clc; close all; 

fname = 'TrajOpt-MovingBody'; 

%% Build Robot 

% file
urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'dev', 'SimpleHAMR', 'urdf', 'HAMRSimple_scaled.urdf');

% options
options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = true;
options.collision = true;

hamr = HAMRSimpleRBM(urdf,options);
v = hamr.constructVisualizer();

%% 
load(fname); 

nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nx = nq+nv;
nu = hamr.getNumInputs();

tt = xtraj.getBreaks();
xx = xtraj.eval(tt);
uu = utraj.eval(tt);

act_dof = hamr.getActuatedJoints();
nq = hamr.getNumPositions();
title_str = {'Front Left Swing', 'Front Left Lift', ...
    'Rear Left Swing', 'Rear Left Lift', ...
    'Front Right Swing', 'Front Right Lift', ...
    'Rear Rear Swing', 'Rear Rear Lift'};

figure(1); clf; hold on;
for i = 1:numel(act_dof)
    subplot(4,2,i); hold on; title(title_str(i))
    yyaxis left; hold on; plot(tt, uu(i,:)); ylabel('Force(N)')
    yyaxis right; hold on; plot(tt, rad2deg(xx(act_dof(i), :))); ylabel('Deflection(deg)')
%     legend('Force', 'Deflection')
end

figure(2); clf;
title_str = {'com x', 'com y', 'com z', 'roll', 'pitch', 'yaw'};
for i = 1:6
    subplot(3,2,i); hold on; title(title_str(i))
    yyaxis left; hold on; plot(tt, xx(i,:)); ylabel('Position')
    yyaxis right; hold on; plot(tt, xx(i+nq, :)); ylabel('Velocity')
%     legend('Position', 'Velocity')
end

phi = zeros(4, numel(tt));
cc = ctraj.eval(ctraj.getBreaks);
for i = 1:numel(tt)
    q = xx(1:nq, i);
    qd = xx(nq+1:nx, i);
    kinsol = hamr.doKinematics(q, qd);
    phi(:,i) = hamr.contactConstraints(kinsol, false);
end


legs = {'FL2', 'RL2', 'FR2', 'RR2'};

figure(3); clf; hold on; 
for i = 1:size(phi, 1)    
    subplot(2,2,i); hold on; title(legs{i}); 
    yyaxis left; plot(tt, phi(i,:)*1e3); ylabel('Distance (um)'); ylim([0, 5])
    yyaxis right; plot(tt, cc(i,:)); ylabel('Force (N)')
%     plot(tt, straj.eval(straj.getBreaks()), 'k')
end


figure(4); clf; hold on; 
for i = 1:size(phi, 1)    
    subplot(2,2,i); hold on; title(legs{i}); 
    plot(tt, phi(i,:).*cc(i,:));
%     yyaxis right; plot(tt, cc(i,:)); ylabel('Force (N)')
    plot(tt, straj.eval(straj.getBreaks()), 'k')
    legend('Constraint', 'S')
end
