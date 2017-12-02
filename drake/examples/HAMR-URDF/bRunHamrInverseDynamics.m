clear; close all; clc;

%% Build Robot

hamr_urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'urdf', 'HAMRSimple_scaled.urdf');

% options
options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = true;       % CHANGE
options.collision = true;

hamr = HAMRSimpleRBM(hamr_urdf,options);

nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nu = hamr.getNumInputs();
nc = hamr.getNumContactPairs();
nd = 4; % pyramidal friction cone approx

%% Load Trajectory
save_dir = '~/Dropbox/CurrentWork/FrictionTrajOpt/MatFiles/TrajOptFiles/';
fname = 'TrajOpt_MovingBody_SimpleSprings7';
td = load([save_dir, fname]);
tt = td.xtraj.getBreaks();
h = mean(diff(tt));

if isempty(td.ctraj)
    td.ctraj = PPTrajectory(zoh(tt, zeros(nc, numel(tt))));
end

if isempty(td.btraj)
    td.btraj = PPTrajectory(zoh(tt, zeros(nc*nd, numel(tt))));
end

xx = td.xtraj.eval(tt);
uu = td.utraj.eval(tt);
cc = td.ctraj.eval(tt);
bb = td.btraj.eval(tt);

%% Form body position

if options.floating
    nfb = 6;
%     xfloating_base = xx(1:nfb, :); 
else
    nfb = 0;
%     xfloating_base = []; %zeros(6, numel(tt)); 
end

%% Form foot positions

leg_links = {'FL2', 'RL2', 'FR2', 'RR2'};     % leg links

pfhamr = [0 0 -14.988382167532292]'; % position of foot in local frame

xfoot_body = zeros(numel(pfhamr), numel(tt), nc);
fkopt = struct();
fkopt.base_or_frame_id = hamr.findLinkId('Chassis'); 
for i = 1:numel(tt)
    kinsol = hamr.doKinematics(xx(1:nq, i), zeros(nv,1));
    for j = 1:nc
        xfoot_body(:,i,j) = hamr.forwardKin(kinsol, hamr.findLinkId(leg_links{j}), pfhamr, fkopt);
    end
end

%% Solve
urdf = {'FL_scaled', 'RL_scaled', 'FR_scaled', 'RR_scaled'};
foot = {'FLL4', 'RLL4', 'FRL4', 'RRL4'};
pfSL = [0, 7.58, -11.35;
    0, 7.58, -11.35;
    0, -7.58, -11.35;
    0, -7.58, -11.35]';

for j = 1:nc
    xfoot_traj = PPTrajectory(foh(tt, xfoot_body(:,:,j)));
    tic;
    [SL.(urdf{j}), xtraj.(urdf{j}),utraj.(urdf{j}),z.(urdf{j}),F.(urdf{j}),info.(urdf{j})]...
        = SLInverseDynDirtran(tt, urdf{j}, foot{j}, pfSL(:,j), xfoot_traj, ...
        td.ctraj(j), td.btraj((j-1)*nd+(1:nd)));
    toc
    
    v = SL.(urdf{j}).constructVisualizer();
    tf = xtraj.(urdf{j}).getBreaks();
    qq = xtraj.(urdf{j}).eval(tf); qq = qq(1:SL.(urdf{j}).getNumPositions(), :);
    qtraj_scaled = PPTrajectory(foh(tf*1e-3, qq));
    qtraj_scaled = qtraj_scaled.setOutputFrame(v.getInputFrame());
    v.playback(qtraj_scaled, struct('slider', true));
    
    % Evaluate final trajectory
    hf = mean(diff(tf));
    xf = xtraj.(urdf{j}).eval(tf);
    vf = xtraj.(urdf{j}).eval(tf+hf/2);
    uf = utraj.(urdf{j}).eval(tf+hf/2);
    cf = utraj.(urdf{j}).eval(tf+hf/2);
    
    % Reconstruct SL foot Traj
    xfootf = zeros(numel(pfSL(:,j)), numel(tf));
    for i = 1:numel(tf)
        qfk = xf(1:SL.(urdf{j}).getNumPositions(), i);
        vfk = xf(SL.(urdf{j}).getNumPositions()+(1:SL.(urdf{j}).getNumVelocities()), i);
        kinsoli = SL.(urdf{j}).doKinematics(qfk, vfk);
        xfootf(:,i) = SL.(urdf{j}).forwardKin(kinsoli, SL.(urdf{j}).findLinkId(foot{j}), pfSL(:,j));
    end
    
    x0 = hamr.getInitialState();
    % Plot
    figure(3*(j-1)+1); clf;
    subplot(2,2,1); hold on; ylabel('Swing Angle')
    plot(tt, rad2deg(xx(nfb+2*j-1,:)));
    plot(tf, rad2deg(xf(3,:)), '--');
    subplot(2,2,2); hold on; ylabel('Lift Angle')
    plot(tt, rad2deg(xx(nfb+2*j,:) - x0(2*j)));
    plot(tf, rad2deg(xf(7,:)), '--');
    subplot(2,2,3); hold on; ylabel('Swing Velocity')
    plot(tt+h/2, rad2deg(xx(nq+nfb+2*j-1,:)));
    plot(tf+h/2, rad2deg(vf(SL.(urdf{j}).getNumPositions()+3,:)), '--');
    subplot(2,2,4); hold on; ylabel('Lift Velocity')
    plot(tt+h/2, rad2deg(xx(nq+nfb+2*j,:)));
    plot(tf+h/2, rad2deg(vf(SL.(urdf{j}).getNumPositions()+7,:)), '--');
    
    figure(3*(j-1)+2); clf;
    subplot(2,2,1); hold on; ylabel('Swing Force')
    plot(tt+h/2, uu(2*j-1,:));
    plot(tf+h/2, uf(1,:), '--');
    subplot(2,2,2); hold on; ylabel('Lift Force')
    plot(tt+h/2, uu(2*j,:));
    plot(tf+h/2, uf(2,:), '--');
    subplot(2,2,3); hold on; ylabel('Normal Force')
    plot(tt+h/2, cc(j,:));
    plot(tf+h/2, cf(SL.(urdf{j}).nu+SL.(urdf{j}).nl+(1:SL.(urdf{j}).nc),:), '--');
    subplot(2,2,4); hold on; ylabel('Tangential Basis')
    plot(tt+h/2, bb(4*(j-1)+(1:4),:), 'b');
    plot(tf+h/2, cf(SL.(urdf{j}).nu+SL.(urdf{j}).nl+SL.(urdf{j}).nc+(1:SL.(urdf{j}).nc*SL.(urdf{j}).nd),:), 'r--');
    
    figure(3*(j-1)+3); clf;
    subplot(3,1,1); hold on; ylabel('Foot X')
    plot(tt, xfoot_body(1,:,j));
    plot(tf, xfootf(1,:), '--');
    
    subplot(3,1,2); hold on; ylabel('Foot Y')
    plot(tt, xfoot_body(2,:,j));
    plot(tf, xfootf(2,:), '--');
    
    subplot(3,1,3); hold on;  ylabel('Foot Z')
    plot(tt, xfoot_body(3,:,j));
    plot(tf, xfootf(3,:), '--');
end
tilefigs;
save([save_dir, fname, '_fullRobot'], 'xtraj', 'utraj')