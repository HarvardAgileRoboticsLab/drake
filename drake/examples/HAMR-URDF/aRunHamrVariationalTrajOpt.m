clear; clc; close all;

load_dir = '~/Dropbox/CurrentWork/FrictionTrajOpt/MatFiles/SimWarmStart/';
save_dir = '~/Dropbox/CurrentWork/FrictionTrajOpt/MatFiles/SimWarmStart/';
fname = 'TROT_0.2N_10Hz';  
sim_traj = load([load_dir, fname, '_Variational2.mat']); 
% sim_tt = sim_traj.tt; 
% sim_hh = mean(diff(sim_tt)); 
% T = sim_tt(end); 
% 
% % last cycle
% tokens = strsplit(fname, '_'); 
% freq = str2double(tokens{3}(1:end-2)); % FOR FIXED BODY
% nc = T*freq*1e-3; 
% nCF0 = find(sim_tt >= (nc-1)/freq/1e-3, 1, 'first');
% nCF1 = numel(sim_tt); 
% 
% % zero out starting x and y pos
% xx_sim = sim_traj.xx(:, nCF0:nCF1);
% xx_sim0 = 0*xx_sim(:,1); xx_sim0(1:2) = xx_sim(1:2, 1); 
% 
% xx_sim = bsxfun(@minus, xx_sim, xx_sim0); 
% 
% traj_init.t = sim_tt(nCF0:nCF1) - sim_tt(nCF0);
% traj_init.x = PPTrajectory(foh(traj_init.t, xx_sim)); 
% traj_init.u = PPTrajectory(zoh(traj_init.t, sim_traj.uu(:, nCF0:nCF1))); 
% traj_init.c = PPTrajectory(zoh(traj_init.t, sim_traj.c_traj(:, nCF0:nCF1))); 
% traj_init.b = PPTrajectory(zoh(traj_init.t, sim_traj.beta_traj(:, nCF0:nCF1))); 
% traj_init.psi = PPTrajectory(zoh(traj_init.t, sim_traj.psi_traj(:, nCF0:nCF1))); 
% traj_init.eta =  PPTrajectory(zoh(traj_init.t, sim_traj.eta_traj(:, nCF0:nCF1))); 
% traj_init.kl =  PPTrajectory(zoh(traj_init.t, sim_traj.kl_traj(:, nCF0:nCF1))); 
traj_init.t = sim_traj.xtraj.getBreaks();
traj_init.x = sim_traj.xtraj;
traj_init.u = sim_traj.utraj;
traj_init.c = sim_traj.ctraj;
traj_init.b = sim_traj.btraj;
traj_init.psi = sim_traj.psitraj;
traj_init.eta = sim_traj.etatraj;
traj_init.kl = sim_traj.kltraj;

%% Optimize

 [hamr,xtraj,utraj,ctraj,btraj,...
    psitraj,etatraj,jltraj, kltraj, straj, ...
    z,F,info,infeasible_constraint_name] = HAMRVariationalPeriodicTrajOptWS(traj_init);

save([save_dir, fname, '_VariationalLHL.mat'], 'xtraj', 'utraj', 'ctraj', 'btraj', 'psitraj', 'etatraj', ...
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
uu = utraj.eval(tt+hh/2);

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

% figure(2); clf;
% title_str = {'com x', 'com y', 'com z', 'roll', 'pitch', 'yaw'};
% for i = 1:6
%     subplot(3,2,i); hold on; title(title_str(i))
%     yyaxis left; hold on; plot(tt, yy(i,:)); ylabel('Position')
%     yyaxis right; hold on; plot(tt+hh/2, vv(i+nq, :)); ylabel('Velocity')
% %     legend('Position', 'Velocity')
% end
% 
% phi = zeros(4, numel(tt));
% cc = ctraj.eval(ctraj.getBreaks + hh/2);
% for i = 2:numel(tt)
%     q = yy(1:nq, i);
%     qd = yy(nq+1:nx, i);
%     kinsol = hamr.doKinematics(q, qd);
%     phi(:,i-1) = hamr.contactConstraints(kinsol, false);
% end


% legs = {'FLL4'};
% 
% figure(3); clf; hold on; 
% for i = 1:size(phi, 1)    
%     subplot(2,2,i); hold on; %title(legs{i}); 
%     yyaxis left; plot(tt, phi(i,:)); ylabel('Distance (mm)'); %ylim([0, 5])
%     yyaxis right; plot(tt, cc(i,:)); ylabel('Force (N)')
% %     plot(tt, straj.eval(straj.getBreaks()), 'k')
% end
% 
% 
% figure(4); clf; hold on; 
% for i = 1:size(phi, 1)    
%     subplot(2,2,i); hold on; %title(legs{i});
%     plot(tt, 1e-4*ones(numel(tt),1), 'k'); 
%     plot(tt, phi(i,:).*cc(i,:)- straj.eval(straj.getBreaks()));
% %     yyaxis right; plot(tt, cc(i,:)); ylabel('Force (N)')
%     
% end
