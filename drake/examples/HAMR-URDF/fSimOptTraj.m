clear; clc; close all;
global u_traj
% kl_traj = [];
u_traj = [];
%% Load Transmission Trajectories

save_dir = '~/Dropbox/CurrentWork/FrictionTrajOpt/MatFiles/TrajOptFiles/';
fname = 'TROT_0.2N_10Hz';
trajTrans = load([save_dir, fname, '_Variational.mat']);

xtrajd = trajTrans.xtraj();
ttd = xtrajd.getBreaks();
hhd = mean(diff(ttd));
xxd = xtrajd.eval(ttd); 


utraj0 = trajTrans.utraj;
uu0 = utraj0.eval(ttd+hhd/2);

% Just inputs
uu = uu0(1:8, :); 
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

v = hamr.constructVisualizer(); 

%% Build OL Inputs and IC

ncyc = 1; 
utraj = PPTrajectory(foh(0:hhd:(ncyc*ttd(end) - hhd), ...
    repmat(uu(:,1:end-1), 1, ncyc)));
utraj = utraj.setOutputFrame(hamr.getInputFrame);

x0 = xtrajd.eval(0);

%% Simulate Open loop
hamr_OL = cascade(utraj, hamr);
xtraj_sim = simulate(hamr_OL, [0 ncyc*ttd(end)], x0);

%% Simulate Closed Loop
kp = 0.5; 
kd = 0.05; %0.3; 

% qa = hamr.getActuatedJoints();
% xtrajA = PPTrajectory(foh(0:hhd:(ncyc*ttd(end) - hhd), ...
%     repmat(xxTransA(:,1:end-1), 1, ncyc)));

PDTracking = HAMRPDTracking(hamr, utraj, xtrajd, kp, kd); 
hamr_CL = feedback(hamr, PDTracking); 

xtraj_sim_CL = simulate(hamr_CL, [0, ncyc*ttd(end)], x0);


%% Plotting

tt_sol = xtraj_sim.getBreaks();
xx_sol_OL = xtraj_sim.eval(tt_sol);
xx_sol_CL = xtraj_sim_CL.eval(tt_sol);

act_dof = hamr.getActuatedJoints();
ndof = hamr.getNumDiscStates();
title_str = {'Front Left Swing', 'Front Left Lift', ...
    'Rear Left Swing', 'Rear Left Lift', ...
    'Front Right Swing', 'Front Right Lift', ...
    'Rear Rear Swing', 'Rear Rear Lift'};

figure(1); clf; hold on;
for i = 1:numel(act_dof)
    subplot(4,2,i); hold on; title(title_str{i})
    plot(tt_sol, u_traj(i, 1:numel(tt_sol)));
    plot(tt_sol, u_traj(i, numel(tt_sol)+(1:numel(tt_sol))));
    plot(ttd+hhd/2, uu(i, :));
    legend('OL Inputs', 'CL Inputs', 'Desired Inputs');
end

% figure(100); clf; hold on;
% for i = 1:numel(act_dof)
%      subplot(4,2,i); hold on; title(title_str{i})
%      yyaxis left; plot(ttd, xxd(act_dof(i), :));
%     yyaxis right; plot(ttd+hhd/2, uu(i, :));
% end

figure(2); clf; hold on;
for i = 1:numel(act_dof)
    subplot(4,2,i); hold on; title(title_str{i})
    plot(tt_sol, xx_sol_OL(act_dof(i), :));
    plot(tt_sol, xx_sol_CL(act_dof(i), :));
    plot(ttd, xxd(act_dof(i), :));
    legend('OL Act Defl', 'CL Act Defl', 'Desired Act Defl');
end

if options.floating
    figure(3); clf;
    title_str = {'Com-X(mm)', 'Com-Y(mm)', 'Com-Z(mm)', 'Roll(deg)', 'Pitch(deg)', 'Yaw(deg)'};
    for i = 1:6
        subplot(2,3,i); hold on; ylabel(title_str(i), 'FontSize', 18)
        xlabel('Time(s)', 'FontSize', 18)
        if i > 3
            plot(tt_sol*1e-3, rad2deg(xx_sol_OL(i,:)), 'LineWidth', 1.5); 
            plot(tt_sol*1e-3, rad2deg(xx_sol_CL(i,:)), 'LineWidth', 1.5); 
            plot(ttd*1e-3, rad2deg(xxd(i,:)), 'LineWidth', 1.5);
            lh = legend('OL', 'CL', 'Desired');
            set(lh, 'box', 'off')
        else
            plot(tt_sol*1e-3, xx_sol_OL(i,:), 'LineWidth', 1.5); 
            plot(tt_sol*1e-3, xx_sol_CL(i,:), 'LineWidth', 1.5); 
            plot(ttd*1e-3, xxd(i,:), 'LineWidth', 1.5);
            lh = legend('OL', 'CL', 'Desired');
            set(lh, 'box', 'off')
        end
    end    
end


pfFull = [0, 7.58, -11.350;
    0, 7.58, -11.350;
    0, -7.58, -11.350;
    0, -7.58, -11.350];

pfCL = zeros([numel(tt_sol), size(pfFull')]);
pfOL = zeros([numel(tt_sol), size(pfFull')]);
pfDesW = zeros([numel(tt_sol), size(pfFull')]);

legs = {'FLL4', 'RLL4', 'FRL4', 'RRL4'};
%
for j = 1:numel(tt_sol)

    qCL = xx_sol_CL(1:ndof/2, j);
    qdCL = xx_sol_CL(ndof/2+1: ndof, j);
    kinsolCL = hamr.doKinematics(qCL, qdCL);
    for i = 1:size(pfFull,1)
        pfCL(j,:,i) = hamr.forwardKin(kinsolCL, hamr.findLinkId(legs{i}), pfFull(i,:)');
    end
    
    qOL = xx_sol_OL(1:ndof/2, j);
    qdOL = xx_sol_OL(ndof/2+1: ndof, j);
    kinsolOL = hamr.doKinematics(qOL, qdOL);
    for i = 1:size(pfFull,1)
        pfOL(j,:,i) = hamr.forwardKin(kinsolOL, hamr.findLinkId(legs{i}), pfFull(i,:)');
    end

    x0 = xtrajd.eval(tt_sol(j));
    q0 = x0(1:ndof/2);
    qd0 =  x0(ndof/2+1:ndof);
    kinsol0 = hamr.doKinematics(q0, qd0);
    for i = 1:size(pfFull,1)
        pfDesW(j,:,i) = hamr.forwardKin(kinsol0, hamr.findLinkId(legs{i}), pfFull(i,:)');
    end
end

figure(4); clf; hold on;
nc = size(pfFull,1);
for i = 1:nc

    subplot(nc,3,3*(i-1)+1); hold on; ylabel('Foot X'); title(legs{i})
    plot(tt_sol, pfOL(:,1,i))
    plot(tt_sol, pfCL(:,1,i))
    plot(tt_sol, pfDesW(:,1,i));
%     legend('Leg Position', 'Desired Leg Position');

    subplot(nc,3,3*(i-1)+2); hold on; ylabel('Foot Y'); title(legs{i})
    plot(tt_sol, pfOL(:,2,i))
    plot(tt_sol, pfCL(:,2,i));
    plot(tt_sol, pfDesW(:,2,i));
%     legend('Leg Position', 'Desired Leg Position');

    subplot(nc,3,3*(i-1)+3); hold on;  ylabel('Foot Z'); title(legs{i})
    plot(tt_sol, pfOL(:,3,i))
    plot(tt_sol, pfCL(:,3,i)); %ylim([0, 1.2*max(pfFullW(:,3,i))]);
    plot(tt_sol, pfDesW(:,3,i));% ylim([0, 1.2*max(pfSimpleW(:,3,i))]);
%     legend('Leg Position', 'Desired Leg Position');

end

tilefigs;

%%
% xtraj_scaled = DTTrajectory(xtraj_sim.getBreaks()*1e-3, xtraj_sim.eval(xtraj_sim.getBreaks()));
% xtraj_scaled = xtraj_scaled.setOutputFrame(xtraj_sim.getOutputFrame());
% options.slider = true;
% v.playback(xtraj_scaled, options);
xtraj_scaled = DTTrajectory(xtraj_sim_CL.getBreaks()*1e-3, xtraj_sim_CL.eval(xtraj_sim_CL.getBreaks()));
xtraj_scaled = xtraj_scaled.setOutputFrame(xtraj_sim_CL.getOutputFrame());
options.slider = true;
v.playback(xtraj_scaled, options);

