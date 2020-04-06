clear; clc; close all;
global lp_b

lp_b = [0, 7.58, -11.350;
    0, 7.58, -11.350;
    0, -7.58, -11.350;
    0, -7.58, -11.350];

%% Load Transmission Trajectories

save_dir = '~/Dropbox/CurrentWork/FrictionTrajOpt/MatFiles/SimWarmStart/';
fname = 'TROT_0.15N_30Hz_TYM_SP';
trajTrans = load([save_dir, fname, '_VariationalSmooth_PlusAct.mat']); %, '_VariationalMU.mat']);
xtrajd = trajTrans.xtraj();
ttd = xtrajd.getBreaks();
hhd = mean(diff(ttd));
xxd = xtrajd.eval(ttd);

% Inputs forces
utrajd = trajTrans.utraj;
uud = utrajd.eval(ttd);

% Input voltages
vtrajd = trajTrans.vtraj; 
vvd = vtrajd.eval(ttd); 
%% Build Full Model

urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2_TYM.urdf');

options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = true; %false;
options.collision = true; %false;
options.dt = 0.4; %0.1; 

% Build robot + visualizer
hamr = HamrTSRBM(urdf, options);
nq = hamr.getNumPositions(); 
nv = hamr.getNumVelocities(); 
nu = hamr.getNumInputs(); 

v = hamr.constructVisualizer(); 

%% Build Actuators
dp.Vb = 225;
dp.Vg = 0;
 
nact = 8;
hr_actuators = HamrActuators(nact, {'FLsact', 'FLlact', 'RLsact', 'RLlact', ...
    'FRsact', 'FRlact', 'RRsact', 'RRlact'}, [], dp);

% for i = 2:2:nact
%     hr_actuators.dummy_bender(i) = hr_actuators.dummy_bender(i).setCFThickness(0.1); 
% end

%% Connect system

%connections from actuators to hamr
hr_actuators = hr_actuators.setOutputFrame(hamr.getInputFrame());
connection1(1).from_output = hr_actuators.getOutputFrame();
connection1(1).to_input = hamr.getInputFrame();

% connections from hamr to actuators
hamr_out = hamr.getOutputFrame();
act_in = hr_actuators.getInputFrame();
act_in = act_in.replaceFrameNum(2, hamr_out.getFrameByName('ActuatorDeflection'));
hr_actuators = hr_actuators.setInputFrame(act_in);
% 
connection2(1).from_output = hamr_out.getFrameByName('ActuatorDeflection');
connection2(1).to_input = act_in.getFrameByName('ActuatorDeflection');

% mimo inputs
input_select(1).system = 1;
input_select(1).input = act_in.getFrameByName('DriveVoltage');

% mimo outputs
output_select(1).system = 2;
output_select(1).output = hamr_out.getFrameByName('HamrPosition');
output_select(2).system = 2;
output_select(2).output = hamr_out.getFrameByName('HamrVelocity');

hamrWact = mimoFeedback(hr_actuators, hamr, connection1, connection2, ...
    input_select, output_select);
%% Build OL Inputs and IC

ncyc = 5; 
ttdN = 0:hhd:(ncyc*ttd(end));

xxdN = [xxd(:,1), repmat(xxd(:,2:end), 1, ncyc)];
if options.floating
   xxCOMdN = reshape(xxdN(1,2:end), [], ncyc);  
   for i = 2:ncyc
    xxCOMdN(:, i) = xxCOMdN(:,i) + xxCOMdN(end, i-1);
   end
   xxdN(1,:) = [0; xxCOMdN(:)];
end

uudN = [uud(:,1), repmat(uud(:,2:end), 1, ncyc)];
utraj = PPTrajectory(zoh(ttdN, uudN)); 
utraj = utraj.setOutputFrame(hamr.getInputFrame);

vvdN = [vvd(:,1), repmat(vvd(:,2:end), 1, ncyc)];
vtraj = PPTrajectory(zoh(ttdN, vvdN)); 
vtraj = vtraj.setOutputFrame(hamrWact.getInputFrame);

x0 = xxd(:,1); 

%% Simulate Open loop

hamr_OL = cascade(vtraj, hamrWact);
xtraj_sim = simulate(hamr_OL, [0 ncyc*ttd(end)], x0);

%% Simulate Closed Loop
% kp = 50; %100; %0.05; %0.5; 
% kd = 40; %20; %0.15; %0.05; %0.3; 
rho = 0.005; 

tracking_opt.ctype = 'actlqr';
tracking_opt.rho = rho; 
% tracking_opt.kp = kp;
% tracking_opt.kd = kd;

% qa = hamr.getActuatedJoints();
% xtrajA = PPTrajectory(foh(0:hhd:(ncyc*ttd(end) - hhd), ...
%     repmat(xxTransA(:,1:end-1), 1, ncyc)));
xtrajd = PPTrajectory(foh(ttdN, xxdN)); 
LegTracking = HAMRLegTracking(hamrWact, hamr, hr_actuators, vtraj, xtrajd, tracking_opt); 

% mimo outputs
output_select(1).system = 1;
output_select(1).output = hamrWact.getOutputFrame.getFrameByName('HamrPosition');
output_select(2).system = 1;
output_select(2).output = hamrWact.getOutputFrame.getFrameByName('HamrVelocity');
output_select(3).system = 2;
output_select(3).output = LegTracking.getOutputFrame();

% hamrWact_CL = mimoFeedback(hamrWact, PDTracking); 
hamr_CL = mimoFeedback(hamrWact, LegTracking, [], [], [], output_select); 
xtraj_sim_CL = simulate(hamr_CL, [0, ncyc*ttd(end)], x0);


%% Plotting

tt_sol = xtraj_sim.getBreaks();
yy_sol_OL = xtraj_sim.eval(tt_sol);
yy_sol_CL = xtraj_sim_CL.eval(tt_sol);

xx_sol_OL = yy_sol_OL(1:nq+nv, :); 
xx_sol_CL = yy_sol_CL(1:nq+nv, :); 

vv_sol_OL = vtraj.eval(tt_sol); 
vv_sol_CL = yy_sol_CL(nq+nv+(1:nu),:); 

act_dof = hamr.getActuatedJoints();
ndof = hamr.getNumDiscStates();
title_str = {'Front Left Swing', 'Front Left Lift', ...
    'Rear Left Swing', 'Rear Left Lift', ...
    'Front Right Swing', 'Front Right Lift', ...
    'Rear Rear Swing', 'Rear Rear Lift'};

figure(1); clf; hold on;
for i = 1:numel(act_dof)
    subplot(4,2,i); hold on; title(title_str{i})
    plot(tt_sol, vv_sol_OL(i,:));
    plot(tt_sol, vv_sol_CL(i,:));
    plot(ttdN+hhd/2, vvdN(i, :));
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
    plot(ttdN, xxdN(act_dof(i), :));
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
            plot(ttdN*1e-3, rad2deg(xxdN(i,:)), 'LineWidth', 1.5);
            lh = legend('OL', 'CL', 'Desired');
            set(lh, 'box', 'off')
        else
            plot(tt_sol*1e-3, xx_sol_OL(i,:), 'LineWidth', 1.5); 
            plot(tt_sol*1e-3, xx_sol_CL(i,:), 'LineWidth', 1.5); 
            plot(ttdN*1e-3, xxdN(i,:), 'LineWidth', 1.5);
            lh = legend('OL', 'CL', 'Desired');
            set(lh, 'box', 'off')
        end
    end    
end

pfCL = zeros([numel(tt_sol), size(lp_b')]);
pfOL = zeros([numel(tt_sol), size(lp_b')]);
pfDesW = zeros([numel(tt_sol), size(lp_b')]);

legs = {'FLL4', 'RLL4', 'FRL4', 'RRL4'};
fkopt.base_or_frame_id = hamr.findLinkId('Chassis');

for j = 1:numel(tt_sol)

    qCL = xx_sol_CL(1:ndof/2, j);
    qdCL = xx_sol_CL(ndof/2+1: ndof, j);
    kinsolCL = hamr.doKinematics(qCL, qdCL);
    for i = 1:size(lp_b,1)
        pfCL(j,:,i) = hamr.forwardKin(kinsolCL, hamr.findLinkId(legs{i}), lp_b(i,:)'); %,fkopt);
    end
    
    qOL = xx_sol_OL(1:ndof/2, j);
    qdOL = xx_sol_OL(ndof/2+1: ndof, j);
    kinsolOL = hamr.doKinematics(qOL, qdOL);
    for i = 1:size(lp_b,1)
        pfOL(j,:,i) = hamr.forwardKin(kinsolOL, hamr.findLinkId(legs{i}), lp_b(i,:)'); %,fkopt);
    end

    xi = xtrajd.eval(tt_sol(j));
    qi = xi(1:ndof/2);
    qdi =  xi(ndof/2+1:ndof);
    kinsol0 = hamr.doKinematics(qi, qdi);
    for i = 1:size(lp_b,1)
        pfDesW(j,:,i) = hamr.forwardKin(kinsol0, hamr.findLinkId(legs{i}), lp_b(i,:)'); %,fkopt);
    end
end

figure(4); clf; hold on;
nc = size(lp_b,1);
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
    legend('OL', 'CL', 'Desired');

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

