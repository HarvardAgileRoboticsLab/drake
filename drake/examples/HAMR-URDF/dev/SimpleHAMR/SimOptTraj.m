clear; clc; close all;

%% Build (open-loop) control input

fname = 'TrajOpt_26-Oct-2017_6';
traj_COM = load(fname);
traj_full = load([fname, '_fullRobot']);

% Build transmission trajectory
xtrajTrans = traj_full.xtraj();
tt = xtrajTrans.FL_scaled.getBreaks();
hh = mean(diff(tt));

xxFL = xtrajTrans.FL_scaled.eval(tt);
xxRL = xtrajTrans.RL_scaled.eval(tt);
xxFR = xtrajTrans.FR_scaled.eval(tt);
xxRR = xtrajTrans.RR_scaled.eval(tt);

nqt = size(xxFL,1)/2;

% Build COM trajectory
xtrajCOM = traj_COM.xtraj();
xxCOM = xtrajCOM.eval(tt);
nqc = size(xxCOM,1)/2;

% Full trajectory
xx = [xxCOM(1:6, :); xxFL(1:nqt, :); xxRL(1:nqt, :); xxFR(1:nqt, :); xxRR(1:nqt, :);
    xxCOM(nqc+(1:6), :); xxFL(nqt+(1:nqt), :); xxRL(nqt+(1:nqt), :); xxFR(nqt+(1:nqt), :); xxRR(nqt+(1:nqt), :)];

% Build input trajectory (included contact and const. forces)
utrajTrans = traj_full.utraj;
uuFull = [utrajTrans.FL_scaled.eval(tt+hh/2);
    utrajTrans.RL_scaled.eval(tt+hh/2);
    utrajTrans.FR_scaled.eval(tt+hh/2);
    utrajTrans.RR_scaled.eval(tt+hh/2)];
nut = size(utrajTrans.FL_scaled.eval(tt), 1);

% Just inputs
uu = [uuFull(1:2, :); uuFull(nut+(1:2),:); uuFull(2*nut+(1:2),:); uuFull(3*nut+(1:2), :)];


%% Model Options

options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = true;
options.collision = true;

%% Build Simpel Model
urdf_simple = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'dev', 'SimpleHAMR', ...
    'urdf', 'HAMRSimple_scaled.urdf');
hamr_simple = HAMRSimpleRBM(urdf_simple, options);
nqS = hamr_simple.getNumPositions();

% Leg trajectory
pfSimple = [0, 0, -14.97;
    0, 0, -14.97;
    0, 0, -14.97;
    0, 0, -14.97];

pfSimpleW = zeros([numel(tt), size(pfSimple')]);

legs = {'FL2', 'RL2', 'FR2', 'RR2'};

for j = 1:numel(tt)
    q = xxCOM(1:nqS, j);
    qd = xxCOM(nqS+(1:nqS), j);
    kinsol = hamr_simple.doKinematics(q, qd);
    for i = 1:size(pfSimple,1)
        pfSimpleW(j,:,i) = hamr_simple.forwardKin(kinsol, hamr_simple.findLinkId(legs{i}), ...
            pfSimple(i,:)');
    end
end


%% Build Full Model

urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2.urdf');

% options
options.z_inactive_guess_tol = 0.1;
options.dt = 0.1; %mean(diff(tt)); %0.1;

% Build robot + visualizer
hamr = HamrTSRBM(urdf, options);
hamr = compile(hamr);
% hamr = hamr.removeAllStateConstraints();
v = hamr.constructVisualizer();

%% Simulate Open loop

ncyc = 1;
tsim = ncyc*tt(end);
if ncyc > 1
    ttsim = linspace(0, tsim, ncyc*(numel(tt)-1));
    utraj = PPTrajectory(foh(ttsim, repmat(uu(:, 1:end-1), 1, ncyc)));
else
    utraj = PPTrajectory(foh(tt, uu));
end
utraj = utraj.setOutputFrame(hamr.getInputFrame);
xtraj = PPTrajectory(foh(tt, xx));
hamr_OL = cascade(utraj, hamr);

x0_hat = xx(:,1);
nq = hamr.getNumPositions();
x0 = resolveConstraints(hamr,x0_hat);
[tf, err_str] = valuecheck(positionConstraints(hamr,double(x0)),zeros(72,1),1e-6);

if tf
    disp('Valid initial condition: simulating...')
    tic;
    xtraj_sim = simulate(hamr_OL, [0 tsim], double(x0));
    tlcp = toc;
    xtraj_scaled = DTTrajectory(xtraj_sim.getBreaks()*1e-3, xtraj_sim.eval(xtraj_sim.getBreaks()));
    xtraj_scaled = xtraj_scaled.setOutputFrame(xtraj_sim.getOutputFrame());
    fprintf('It took %fs to simulate %fs of realtime. \nThats %fx \n', ...
        tlcp, tsim/1000, 1000*tlcp/tsim)
    options.slider = true;
    %     xtraj.tt = xtraj.tt/1000;
    v.playback(xtraj_scaled, options);
else
    disp('invalid initial condition...')
end

%% Plotting

tt_sol = xtraj_sim.getBreaks();
xx_sol = xtraj_sim.eval(tt_sol);
uu_sol = utraj.eval(tt_sol);

act_dof = hamr.getActuatedJoints();
ndof = hamr.getNumDiscStates();
title_str = {'Front Left Swing', 'Front Left Lift', ...
    'Rear Left Swing', 'Rear Left Lift', ...
    'Front Right Swing', 'Front Right Lift', ...
    'Rear Rear Swing', 'Rear Rear Lift'};

figure(2); clf; hold on;
for i = 1:numel(act_dof)
    subplot(4,2,i); hold on; title(title_str(i))
    yyaxis left; hold on; plot(tt_sol, uu_sol(i,:), 'b')
    yyaxis right; hold on; plot(tt_sol, xx_sol(act_dof(i), :)*1e3)
    legend('Force', 'Deflection')
end

pfFull = [0, 7.58, -11.350;
    0, 7.58, -11.350;
    0, -7.58, -11.350;
    0, -7.58, -11.350];

pfFullW = zeros([numel(tt_sol), size(pfFull')]);

legs = {'FLL4', 'RLL4', 'FRL4', 'RRL4'};

for j = 1:numel(tt_sol)
    q = xx_sol(1:ndof/2, j);
    qd = xx_sol(ndof/2+1: ndof, j);
    kinsol = hamr.doKinematics(q, qd);
    for i = 1:size(pfFull,1)
        pfFullW(j,:,i) = hamr.forwardKin(kinsol, hamr.findLinkId(legs{i}), pfFull(i,:)');
    end
end

figure(3); clf; hold on;

nc = size(pfFull,1);
for i = 1:nc
    
    subplot(nc,3,3*(i-1)+1); hold on; ylabel('Foot X'); title(legs{i})
    plot(tt_sol, pfFullW(:,1,i))
    plot(tt, pfSimpleW(:,1,i));
    
    subplot(nc,3,3*(i-1)+2); hold on; ylabel('Foot Y'); title(legs{i})
    plot(tt_sol, pfFullW(:,2,i))
    plot(tt, pfSimpleW(:,2,i));
    
    subplot(nc,3,3*(i-1)+3); hold on;  ylabel('Foot Z'); title(legs{i})
    yyaxis left;  plot(tt_sol, pfFullW(:,3,i)); ylim([0, 1.2*max(pfFullW(:,3,i))]);
    yyaxis right; plot(tt, pfSimpleW(:,3,i)); ylim([0, 1.2*max(pfSimpleW(:,3,i))]);
    
end


figure(4); clf;
title_str = {'com x', 'com y', 'com z', 'roll', 'pitch', 'yaw'};
for i = 1:6
    subplot(2,3,i); hold on; title(title_str(i))
    if i > 3
        plot(tt_sol, rad2deg(xx_sol(i,:))); plot(tt, rad2deg(xxCOM(i,:)));
    else
        plot(tt_sol, xx_sol(i,:)); plot(tt, xxCOM(i,:));
    end
    
    %     yyaxis right; hold on; plot(tt_sol, xx_sol(i+ndof/2, :)); plot(tt, xxCOM(nqc+i, :));
end

contact_opt.use_bullet = false;
phi = zeros(4, numel(tt_sol));
for i = 1:numel(tt_sol)
    q = xx_sol(1:ndof/2, i);
    qd = xx_sol(ndof/2+1:ndof, i);
    kinsol = hamr.doKinematics(q, qd);
    phi(:,i) = hamr.contactConstraints(kinsol, false, contact_opt);
end

%     figure(5); clf;
%     for i = 1:size(phi, 1)
%         subplot(2,2,i); hold on;
%         yyaxis right; plot(tt, phi(i,:)); ylabel('Leg Height')
%         yyaxis left; plot(tt, c_traj(i,:)); ylabel('Force')
%         plot(tt, beta_traj(4*(i-1)+(1:4), :), 'k');
%     end

    
% %% saving
% savedir = '';
%
% if SAVE_FLAG
%     fname = [gait, '_', num2str(dp.Vb) 'V_', num2str(1e3*fd), 'Hz_', num2str(options.mu*100), '.mat'];
%     disp(['Saving as ', fname]);
%     save([savedir, fname], 'tt', 'xx', 'uu', 'kl_traj', 'jl_traj', 'c_traj', 'beta_traj', 'eta_traj', 'psi_traj', 'Vact');
% end
