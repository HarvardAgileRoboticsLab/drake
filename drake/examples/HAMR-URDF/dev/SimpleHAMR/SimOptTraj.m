clear; clc; close all;
global u_traj kl_traj

kl_traj = [];
u_traj = [];
%% Load Transmission Trajectories

fname = 'TrajOpt_MovingBody';
trajTrans = load([fname, '_fullRobot']);
xtrajTrans = trajTrans.xtraj();
ttTrans = xtrajTrans.FL_scaled.getBreaks();
hhTrans = mean(diff(ttTrans));

% Build transmission trajectory
xtrajTrans = [xtrajTrans.FL_scaled(1:8);
    xtrajTrans.RL_scaled(1:8);
    xtrajTrans.FR_scaled(1:8);
    xtrajTrans.RR_scaled(1:8); 
    xtrajTrans.FL_scaled(8+(1:8));
    xtrajTrans.RL_scaled(8+(1:8));
    xtrajTrans.FR_scaled(8+(1:8));
    xtrajTrans.RR_scaled(8+(1:8))];
xxTrans = xtrajTrans.eval(ttTrans);

% Build input trajectory (included contact and const. forces)
utrajTrans = trajTrans.utraj;
uuTrans = [utrajTrans.FL_scaled.eval(ttTrans+hhTrans/2);
    utrajTrans.RL_scaled.eval(ttTrans+hhTrans/2);
    utrajTrans.FR_scaled.eval(ttTrans+hhTrans/2);
    utrajTrans.RR_scaled.eval(ttTrans+hhTrans/2)];

% Just inputs
nut = size(utrajTrans.FL_scaled.eval(ttTrans), 1);
uu = [uuTrans(1:2, :); uuTrans(nut+(1:2),:); uuTrans(2*nut+(1:2),:); uuTrans(3*nut+(1:2), :)];

%% Build Full Model

urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2.urdf');

options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = true; %false;
options.collision = true; %false;
options.dt = 0.5; %hhTrans; %0.1; 

% Build robot + visualizer
hamr = HamrTSRBM(urdf, options);
hamr = compile(hamr);
v = hamr.constructVisualizer();

% % remove joint limits
% [lb, ub] = hamr.getJointLimits(); 
% hamr = hamr.setJointLimits(-Inf*lb, Inf*ub); 
% hamr = compile(hamr); 

%% Load COM Trajectories
if options.floating
    trajCOM = load(fname);
    ttCOM = trajCOM.xtraj.getBreaks();
    xtrajCOM = trajCOM.xtraj;
    xxCOM = xtrajCOM.eval(ttCOM);
    xxCOM = xxCOM([1:6, 14+(1:6)], :);
end

%% Build OL Inputs and IC

ncyc = 1; 
utraj = PPTrajectory(zoh(0:hhTrans:(ncyc*ttTrans(end) - hhTrans), ...
    repmat(uu(:,1:end-1), 1, ncyc)));
utraj = utraj.setOutputFrame(hamr.getInputFrame);

x0 = xtrajTrans.eval(0);
if options.floating
    x0COM = xtrajCOM.eval(0);
    x0 = [x0COM(1:6); x0(1:32); x0COM(14+(1:6));x0(32+(1:32))];
end

%% Simulate Open loop
hamr_OL = cascade(utraj, hamr);
xtraj_sim = simulate(hamr_OL, [0 ncyc*ttTrans(end)], x0);

%% Simulate Closed Loop
kp = 10;
kd = 1; 

hamr_CL = hamr.pdtrajtracking(kp, kd);

xxTransA = xxTrans(hamr.getActuatedJoints()-6, :);
xtrajA = PPTrajectory(foh(0:hhTrans:(ncyc*ttTrans(end) - hhTrans), ...
    repmat(xxTransA(:,1:end-1), 1, ncyc)));
xtraj_sim_CL = simulate(cascade(setOutputFrame(xtrajA, hamr_CL.getInputFrame), hamr_CL), ...
    [0, ncyc*ttTrans(end)], x0);


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
    plot(ttTrans+hhTrans/2, uu(i, :));
    legend('OL Inputs', 'CL Inputs', 'Desired Inputs');
end


figure(2); clf; hold on;
for i = 1:numel(act_dof)
    subplot(4,2,i); hold on; title(title_str{i})
    plot(tt_sol, xx_sol_OL(act_dof(i), :));
    plot(tt_sol, xx_sol_CL(act_dof(i), :));
    plot(ttTrans, xxTrans(act_dof(i)-6, :));
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
            plot(ttCOM*1e-3, rad2deg(xxCOM(i,:)), 'LineWidth', 1.5);
            lh = legend('OL', 'CL', 'Desired');
            set(lh, 'box', 'off')
        else
            plot(tt_sol*1e-3, xx_sol_OL(i,:), 'LineWidth', 1.5); 
            plot(tt_sol*1e-3, xx_sol_CL(i,:), 'LineWidth', 1.5); 
            plot(ttCOM*1e-3, xxCOM(i,:), 'LineWidth', 1.5);
            lh = legend('OL', 'CL', 'Desired');
            set(lh, 'box', 'off')
        end
    end    
end

tilefigs;

%%
xtraj_scaled = DTTrajectory(xtraj_sim_CL.getBreaks()*1e-3, xtraj_sim_CL.eval(xtraj_sim_CL.getBreaks()));
xtraj_scaled = xtraj_scaled.setOutputFrame(xtraj_sim_CL.getOutputFrame());
options.slider = true;
v.playback(xtraj_scaled, options);

%%


% pfFull = [0, 7.58, -11.350;
%     0, 7.58, -11.350;
%     0, -7.58, -11.350;
%     0, -7.58, -11.350];
%
% pfSimW = zeros([numel(tt_sol), size(pfFull')]);
% pfDesW = zeros([numel(tt_sol), size(pfFull')]);
%
% legs = {'FLL4', 'RLL4', 'FRL4', 'RRL4'};
% %
% for j = 1:numel(tt_sol)
%
%     q = xx_sol(1:ndof/2, j);
%     qd = xx_sol(ndof/2+1: ndof, j);
%     kinsol = hamr.doKinematics(q, qd);
%     for i = 1:size(pfFull,1)
%         pfSimW(j,:,i) = hamr.forwardKin(kinsol, hamr.findLinkId(legs{i}), pfFull(i,:)');
%     end
%
%     xdTrans = xtrajTrans.eval(tt_sol(j));
%     xdCOM = xtrajCOM.eval(tt_sol(j));
%     q0 = [xdCOM(1:6); xdTrans(1:32)];
%     qd0 = [xdCOM(14+(1:6)); xdTrans(32+(1:32))];
%     kinsol0 = hamr.doKinematics(q0, qd0);
%     for i = 1:size(pfFull,1)
%         pfDesW(j,:,i) = hamr.forwardKin(kinsol0, hamr.findLinkId(legs{i}), pfFull(i,:)');
%     end
% end
%
% figure(3); clf; hold on;
% nc = size(pfFull,1);
% for i = 1:nc
%
%     subplot(nc,3,3*(i-1)+1); hold on; ylabel('Foot X'); title(legs{i})
%     plot(tt_sol, pfSimW(:,1,i))
%     plot(tt_sol, pfDesW(:,1,i), '--');
%     legend('Leg Position', 'Desired Leg Position');
%
%     subplot(nc,3,3*(i-1)+2); hold on; ylabel('Foot Y'); title(legs{i})
%     plot(tt_sol, pfSimW(:,2,i));
%     plot(tt_sol, pfDesW(:,2,i), '--');
%     legend('Leg Position', 'Desired Leg Position');
%
%     subplot(nc,3,3*(i-1)+3); hold on;  ylabel('Foot Z'); title(legs{i})
%     plot(tt_sol, pfSimW(:,3,i)); %ylim([0, 1.2*max(pfFullW(:,3,i))]);
%     plot(tt_sol, pfDesW(:,3,i), '--');% ylim([0, 1.2*max(pfSimpleW(:,3,i))]);
%     legend('Leg Position', 'Desired Leg Position');
%
% end
