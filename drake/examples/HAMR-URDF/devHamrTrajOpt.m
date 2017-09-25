clear; clc; close all;

LOAD = 1; 
[hamr,xtraj,utraj,ctraj,btraj,...
    psitraj,etatraj,jltraj, kltraj, straj,z,F,info,infeasible_constraint_name]= runVariationalTrajOpt(1);

save(['TrajOpt_', datestr(datetime)], 'hamr', 'xtraj', 'utraj', 'ctraj', 'btraj', 'psitraj', 'etatraj', ...
    'jltraj', 'kltraj', 'straj'); 

%% Playback

v = hamr.constructVisualizer();
qq = xtraj.eval(xtraj.getBreaks()); qq = qq(1:50, :); 
qtraj_scaled = PPTrajectory(foh(xtraj.getBreaks()*1e-3, qq));
qtraj_scaled = qtraj_scaled.setOutputFrame(v.getInputFrame());
v.playback(qtraj_scaled, struct('slider', true));


%% Plotting
% load('NomTraj.mat')
% ttnom = tt; 
% uunom = yy(101:108,:);
% 
% urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMRVariational_scaledV2_LB.urdf');
% 
% % options
% options.terrain = RigidBodyFlatTerrain();
% options.ignore_self_collisions = true;
% options.collision_meshes = false;
% options.use_bullet = false;
% options.floating = true;
% options.collision = true;


% hamr = HamrVariationalRBM(urdf,options);

% load TrajOpt_21-Sep-2017 17:16:44

nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nx = nq+nv;
nu = hamr.getNumInputs();

tt = xtraj.getBreaks();
yy = xtraj.eval(tt);
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
    yyaxis right; hold on; plot(tt, yy(act_dof(i), :)); ylabel('Deflection(mm)')
%     legend('Force', 'Deflection')
end

figure(2); clf;
title_str = {'com x', 'com y', 'com z', 'roll', 'pitch', 'yaw'};
for i = 1:6
    subplot(3,2,i); hold on; title(title_str(i))
    yyaxis left; hold on; plot(tt, yy(i,:)); ylabel('Position')
    yyaxis right; hold on; plot(tt, yy(i+nq, :)); ylabel('Velocity')
%     legend('Position', 'Velocity')
end

phi = zeros(4, numel(tt));
cc = ctraj.eval(ctraj.getBreaks);
for i = 1:numel(tt)
    q = yy(1:nq, i);
    qd = yy(nq+1:nx, i);
    kinsol = hamr.doKinematics(q, qd);
    phi(:,i) = hamr.contactConstraints(kinsol, false);
end


legs = {'FLL4', 'RLL4', 'FRL4', 'RRL4'};

figure(3); clf; hold on; 
for i = 1:size(phi, 1)    
    subplot(2,2,i); hold on; title(legs{i}); 
    yyaxis left; plot(tt, phi(i,:)); ylabel('Distance (mm)'); ylim([0, 5])
    yyaxis right; plot(tt, cc(i,:)); ylabel('Force (N)')
    plot(tt, straj.eval(straj.getBreaks()), 'k')
end


figure(4); clf; hold on; 
for i = 1:size(phi, 1)    
    subplot(2,2,i); hold on; title(legs{i}); 
    plot(tt, phi(i,:).*cc(i,:));
%     yyaxis right; plot(tt, cc(i,:)); ylabel('Force (N)')
    plot(tt, straj.eval(straj.getBreaks()), 'k')
end
% lp_b = [0, 7.540, -11.350;
%     0, 7.540, -11.350;
%     0, -7.540, -11.350;
%     0, -7.540, -11.350];
% 
% 
% lp_g = zeros([numel(tt), size(lp_b')]);
% 
% for j = 1:numel(tt)
%     q = yy(1:nq, i);
%     qd = yy(nq+1:nx, i);
%     kinsol = hamr.doKinematics(q, qd);
%     for i = 1:size(lp_b,1)
%         lp_g(j,:,i) = hamr.forwardKin(kinsol, hamr.findLinkId(legs{i}),lp_b(i,:)');
%     end
% end
% figure(4); clf; hold on;
% title('Leg XZ')
% for i = 1:size(lp_b,1)    
%     plot(tt, lp_g(:,3,i))
% end
% legend(legs)
% axis equal; 
