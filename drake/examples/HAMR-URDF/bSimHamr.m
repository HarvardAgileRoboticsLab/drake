clear; clc; close all;
global kl_traj jl_traj c_traj beta_traj psi_traj eta_traj
save_dir = '~/Dropbox/CurrentWork/FrictionTrajOpt/MatFiles/AgileBehaviors/';

%% Load Rigid Body

urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2_TYM.urdf');

kl_traj = [];
jl_traj = [];
c_traj = [];
beta_traj = [];
psi_traj = [];
eta_traj = [];

% options
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.z_inactive_guess_tol = 0.1;
options.use_bullet = false;

% options to change
options.dt = 1;
gait = 'TROT';
SAVE_FLAG = 0;
ISFLOAT = true; % floating (gnd contact) or in air (not floating)

if ISFLOAT
    options.floating = ISFLOAT;
    options.collision = ISFLOAT;
    x0 = zeros(76, 1); x0(3) = 12.69;
    options.terrain = RigidBodyFlatTerrain();
    
else
    options.floating = ISFLOAT;
    options.collision = ISFLOAT;
    x0 = zeros(64, 1);
    options.terrain = [];
end

% Build robot + visualizer
hamr = HamrTSRBM(urdf, options);
hamr = compile(hamr);
v = hamr.constructVisualizer();

%% Build (open-loop) control input

fd = 0.01;         % drive frequency (Hz)
tsim = 500;
FmaxS = 0.25;
FmaxL = 0.35;


t = 0:options.dt:tsim;

switch gait
    case 'SINGLE_JUMP'
        t0 = 20; ind0 = find(t >= t0, 1, 'first');
        t1 = 35; ind1 = find(t > t1, 1, 'first');
        t2 = 45; ind2 = find(t > t2, 1, 'first');
%         t3 = 100; ind3 = find(t > t3, 1, 'first');
%         t4 = 120; ind4 = find(t > t4, 1, 'first');
        
        %         t3 = 65; ind3 = find(t > t3, 1, 'first');
        
        uu = zeros(8, numel(t));
        uu(2:2:end, 1:ind0) = repmat((-0.8*FmaxL/t0)*t(1:ind0), 4, 1);         % pull lifts up
        uu(2:2:end, ind0+1:ind1) = -0.8*FmaxL;                                 % hold
        uu(2:2:end, ind1+1:ind2) = repmat(1.8*FmaxL/(t2-t1)*(t(ind1+1:ind2) - t(ind1)) ...
            - 0.8*FmaxL, 4, 1);         % bring down rapidly
        uu(2:2:end, ind2+1:end) = FmaxL;            % hold at bottom
%         uu(2:2:end, ind3+1:ind4) = repmat(-FmaxL/(t4-t3)*(t(ind3+1:ind4) - t(ind3)) ...
%             + FmaxL, 4, 1);                % hold at bottom

        uu([2, 4], :) = -uu([2, 4], :);             % flip signs for left side
        
    case 'WHEELI'
        t0 = 20; ind0 = find(t >= t0, 1, 'first');
        t0a = 30; ind0a = find(t >= t0a, 1, 'first');
        t1 = 35; ind1 = find(t > t1, 1, 'first');
        t1a = 35; ind1a = find(t >= t1a, 1, 'first');
        t2 = 55; ind2 = find(t > t2, 1, 'first');
        t2a = 55; ind2a = find(t >= t2a, 1, 'first');
        
        twl = 100; indw = find(t >= twl, 1, 'first');
        tws = twl + 0.25/fd; indws = find(t >= tws, 1, 'first');
        
        
        uu = zeros(8, numel(t));
        
        uu([2,6], 1:ind0) = repmat((-FmaxL/t0)*t(1:ind0), 2, 1);         % pull front lifts up
        uu([1,5],1:ind0a) = repmat((FmaxS/t0a)*t(1:ind0a), 2, 1);      % pull swings back
        
        uu([2,6], ind0+1:ind1) = -FmaxL;                                 % hold
        uu([1,5],ind0a+1:end) = FmaxS;
        
        uu([2,6], ind1+1:ind2) = repmat(2*FmaxL/(t2-t1)*(t(ind1+1:ind2) - t(ind1)) - FmaxL, 2, 1);         % bring down rapidly        
        uu([1,5], ind1a+1:ind2a) = repmat(-2*FmaxS/(t2a-t1a)*(t(ind1a+1:ind2a) - t(ind1a)) + FmaxS, 2, 1);         % bring forward rapidly
%         uu([3,7], ind1+1:ind2) = repmat(FmaxS/(t2-t1)*(t(ind1+1:ind2) - t(ind1)), 2, 1);  
        
        uu([2,6], ind2+1:end) = FmaxL;            % hold at bottom
        uu([1,5], ind2a+1:end) = -FmaxS;            % hold at front
  
        uu([2, 4], :) = -uu([2, 4], :);             % flip signs for left side
        uu([3, 5], :) = -uu([3, 5], :);             % flip signs for left side
        
        uu(4, indw+1:end) = 0.3*FmaxL*sin(2*pi*fd*(t(indw+1:end)-twl)); 
        uu(8, indw+1:end) = 0.3*FmaxL*sin(2*pi*fd*(t(indw+1:end)-twl) + pi);
        
        uu(3, indws+1:end) = FmaxS*sin(2*pi*fd*(t(indws+1:end)-tws)); 
        uu(7, indws+1:end) = FmaxS*sin(2*pi*fd*(t(indws+1:end)-tws) + pi);
        
%         uu(4, :) = FmaxL*sin(2*pi*fd*(t-t2) + pi); 
%         uu(8, :) = FmaxL*sin(2*pi*fd*(t-t2) + pi);
        
        
    case 'TROT'
        uu = [FmaxS*sin(2*pi*fd*t + pi/2);            % FLswing
            FmaxL*sin(2*pi*fd*t);                       % FLlift
            FmaxS*sin(2*pi*fd*t + pi/2);              % RLSwing
            FmaxL*sin(2*pi*fd*t + pi);                       % RLLift
            FmaxS*sin(2*pi*fd*t + pi/2);                % FRswing
            FmaxL*sin(2*pi*fd*t);                       % FRlift
            FmaxS*sin(2*pi*fd*t + pi/2);                % RRSwing
            FmaxL*sin(2*pi*fd*t + pi)];               % RRLift
    case 'PRONK'
        uu = [FmaxS*sin(2*pi*fd*t + pi/2);            % FLswing
            FmaxL*sin(2*pi*fd*t);                       % FLlift
            FmaxS*sin(2*pi*fd*t - pi/2);              % RLSwing
            FmaxL*sin(2*pi*fd*t);                       % RLLift
            FmaxS*sin(2*pi*fd*t - pi/2);              % FRswing
            FmaxL*sin(2*pi*fd*t + pi);                  % FRlift
            FmaxS*sin(2*pi*fd*t + pi/2);              % RRSwing
            FmaxL*sin(2*pi*fd*t + pi)];                      % RRLift
    otherwise
        uu = zeros(8, numel(t));
end

% ramp
tramp = 2/fd;
ramp = t/tramp; ramp(t >= tramp) = 1;

uu = bsxfun(@times, ramp, uu); %+ 0.5*(dp.Vb - dp.Vg);
u = PPTrajectory(foh(t, uu));
u = setOutputFrame(u, hamr.getInputFrame());

figure(1); clf;
plot(t, uu(1,:), t, uu(2,:), '--');
plot(t, uu(3,:), t, uu(4,:), '--');

legend('Swing Drive', 'Lift Drive')

%% Simulate Open loop

hamr_OL = cascade(u, hamr);
nQ = hamr.getNumPositions();
x0_hat = hamr.getManipulator.positionConstraints(x0);
[tf, err_str] = valuecheck(positionConstraints(hamr,x0),zeros(72,1),1e-6);

% x0_hat = hamr.positionConstraints(x0(1:nQ));
% [tf, err_str] = valuecheck(positionConstraints(hamr,x0(1:nQ)),zeros(72,1),1e-6);

% tf = 1;
if tf
    disp('Valid initial condition: simulating...')
    tic;
    xtraj = simulate(hamr_OL, [0 tsim], x0);
    tlcp = toc;
    xtraj_scaled = DTTrajectory(xtraj.getBreaks()*1e-3, xtraj.eval(xtraj.getBreaks()));
    xtraj_scaled = xtraj_scaled.setOutputFrame(xtraj.getOutputFrame());
    fprintf('It took %fs to simulate %fs of realtime. \nThats %fx \n', ...
        tlcp, tsim/1000, 1000*tlcp/tsim)
    options.slider = true;
    %     xtraj.tt = xtraj.tt/1000;
    v.playback(xtraj_scaled, options);
else
    disp('invalid initial condition...')
end

%% Plotting
tt = xtraj.getBreaks();
yy = xtraj.eval(tt);
xx = yy(1:2*nQ,:);
% uu = yy(2*nQ+1:end,:);

act_dof = hamr.getActuatedJoints();
ndof = hamr.getNumDiscStates();
title_str = {'Front Left Swing', 'Front Left Lift', ...
    'Rear Left Swing', 'Rear Left Lift', ...
    'Front Right Swing', 'Front Right Lift', ...
    'Rear Rear Swing', 'Rear Rear Lift'};

figure(2); clf; hold on;
for i = 1:numel(act_dof)
    subplot(4,2,i); hold on; title(title_str(i))
    yyaxis left; hold on; plot(tt, uu(i,:))
    %     yyaxis left; plot(tt, yy(act_dof(i), :)*1e3);
    %     yyaxis right; plot(tt, Vact(i,:));
    %     legend('Deflection', 'Force')
    yyaxis right; hold on; plot(tt, yy(act_dof(i), :)*1e3)
    legend('Force', 'Deflection')
end

lp_b = [0, 7.540, -11.350;
    0, 7.540, -11.350;
    0, -7.540, -11.350;
    0, -7.540, -11.350];

lp_g = zeros([numel(t), size(lp_b')]);

legs = {'FLL4', 'RLL4', 'FRL4', 'RRL4'};

for j = 1:numel(tt)
    q = yy(1:ndof/2, j);
    qd = yy(ndof/2+1: ndof, j);
    kinsol = hamr.doKinematics(q, qd);
    for i = 1:size(lp_b,1)
        lp_g(j,:,i) = hamr.forwardKin(kinsol, hamr.findLinkId(legs{i}), lp_b(i,:)');
    end
end

figure(3); clf; hold on;
for i = 1:size(lp_b,1)
    %     subplot(2,2,i); hold on; title(legs{i});
    plot((lp_g(:,1,i) - mean(lp_g(:,1,i))), ...
        (lp_g(:,3,i) - mean(lp_g(:,3,i))))
    axis equal;
    %     axis([-2.5, 2.5, -2.5, 2.5])
end
legend(legs)

if ISFLOAT
    figure(4); clf;
    title_str = {'com x', 'com y', 'com z', 'roll', 'pitch', 'yaw'};
    for i = 1:6
        subplot(3,2,i); hold on; title(title_str(i))
        yyaxis left; hold on; plot(tt, yy(i,:))
        yyaxis right; hold on; plot(tt, yy(i+ndof/2, :))
    end
    
    contact_opt.use_bullet = false;
    phi = zeros(4, numel(tt));
    for i = 1:numel(tt)
        q = yy(1:ndof/2, i);
        qd = yy(ndof/2+1:ndof, i);
        kinsol = hamr.doKinematics(q, qd);
        phi(:,i) = hamr.contactConstraints(kinsol, false, contact_opt);
    end
    
    figure(5); clf;
    for i = 1:size(phi, 1)
        subplot(2,2,i); hold on;
        yyaxis right; plot(tt, phi(i,:)); ylabel('Leg Height')
        yyaxis left; plot(tt, c_traj(i,:)); ylabel('Force')
        plot(tt, beta_traj(4*(i-1)+(1:4), :), 'k');
    end
end

%% saving

if SAVE_FLAG
    fname = [gait, '_', num2str(FmaxS), 'N_', num2str(1e3*fd), 'Hz'];
    disp(['Saving as ', fname]);
    save([save_dir, fname, '_TYM.mat'], 'tt', 'xx', 'uu', 'kl_traj', 'jl_traj', 'c_traj', 'beta_traj', 'eta_traj', 'psi_traj');
end
