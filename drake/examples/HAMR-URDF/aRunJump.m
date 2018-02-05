clear; clc; close all; 

load_dir = '~/Dropbox/CurrentWork/FrictionTrajOpt/MatFiles/AgileBehaviors/';
save_dir = '~/Dropbox/CurrentWork/FrictionTrajOpt/MatFiles/AgileBehaviors/';
fname = 'SINGLE_JUMP_0.25N_10Hz_TYM';

sim_traj = load([load_dir, fname, '.mat']);

traj_init.t = sim_traj.tt; 
traj_init.x = PPTrajectory(foh(traj_init.t, sim_traj.xx)); 
traj_init.u = PPTrajectory(zoh(traj_init.t, sim_traj.uu)); 
traj_init.c = PPTrajectory(zoh(traj_init.t, sim_traj.c_traj)); 
traj_init.b = PPTrajectory(zoh(traj_init.t, sim_traj.beta_traj)); 
traj_init.psi = PPTrajectory(zoh(traj_init.t, sim_traj.psi_traj)); 
traj_init.eta =  PPTrajectory(zoh(traj_init.t, sim_traj.eta_traj)); 
traj_init.kl =  PPTrajectory(zoh(traj_init.t, sim_traj.kl_traj)); 


[hamr,xtraj,utraj,ctraj,btraj,...
    psitraj,etatraj,jltraj, kltraj, straj, ...
    z,F,info,infeasible_constraint_name] = HAMRVariationalStepUp(traj_init); 

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
vv = xtraj.eval(tt);
uu = utraj.eval(tt);

act_dof = hamr.getActuatedJoints();
nq = hamr.getNumPositions();

figure(1); clf; hold on;
for i = 1:nu
    subplot(4,2,i); hold on; %title(title_str(i))
    yyaxis left; hold on; plot(tt, uu(i,:)); ylabel('Force(N)')
    yyaxis right; hold on; plot(tt, yy(act_dof(i), :)); ylabel('Deflection(mm)')
    yyaxis right; hold on; plot(tt+hh/2, vv(nq+act_dof(i), :));
    legend('Force', 'Deflection', 'Deflection Rate')
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

%%
save([save_dir, fname, '4.mat'], 'xtraj', 'utraj', ...
    'ctraj', 'btraj', 'psitraj', 'etatraj',  'jltraj', 'kltraj', 'straj');


