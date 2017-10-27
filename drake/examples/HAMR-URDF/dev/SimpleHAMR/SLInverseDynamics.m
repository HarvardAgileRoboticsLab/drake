clear; close all; clc; 

%% Build Robot

hamr_urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'dev', 'SimpleHAMR', ...
    'urdf', 'HAMRSimple_scaled.urdf');

% options
options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = true;
options.collision = true;

hamr = HAMRSimpleRBM(hamr_urdf,options);

nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nc = hamr.getNumContactPairs();
nd = 4; % pyramidal friction cone approx

%% Load Trajectory

fname = 'TrajOpt-MovingBody';
td = load(fname);
tt = td.xtraj.getBreaks(); 
h = mean(diff(tt));

if isempty(td.ctraj)
    td.ctraj = PPTrajectory(zoh(tt, zeros(nc, numel(tt)))); 
end

if isempty(td.btraj)
    td.btraj = PPTrajectory(zoh(tt, zeros(nc*nd, numel(tt)))); 
end

xx = td.xtraj.eval(tt);
uu = td.utraj.eval(tt+h/2);
cc = td.ctraj.eval(tt);
bb = td.btraj.eval(tt);

%% Form foot position and Forces

pfhamr = [0 0 -14.97]'; % position of foot in local frame
xfoot_body = zeros(numel(pfhamr), numel(tt));

fkopt.base_or_frame_id = hamr.findLinkId('Chassis');

for i = 1:numel(tt)
    kinsol = hamr.doKinematics(xx(1:nq, i), zeros(nv,1));
    xfoot_body(:,i) = hamr.forwardKin(kinsol, hamr.findLinkId('FL2'), pfhamr, fkopt);
end

%% Solve
xfoot_traj = PPTrajectory(foh(tt, xfoot_body));
tic;
pfSL = [0; 7.58; -11.35];
urdfSL = 'FL_scaled'; 
footSL = 'FLL4';
[SL, xtraj,utraj,z,F,info] = SLInverseDynDirtran(tt, urdfSL, footSL, pfSL, xfoot_traj, ...
    td.ctraj(1), td.btraj(1:nd));
toc

%% Animate and Plot
v = SL.constructVisualizer(); 
tf = xtraj.getBreaks();  
qq = xtraj.eval(tf); qq = qq(1:SL.getNumPositions(), :); 
qtraj_scaled = PPTrajectory(foh(tf*1e-3, qq));
qtraj_scaled = qtraj_scaled.setOutputFrame(v.getInputFrame());
v.playback(qtraj_scaled, struct('slider', true));

% Evaluate final trajectory
hf = mean(diff(tf)); 
xf = xtraj.eval(tf);
vf = xtraj.eval(tf+hf/2);
uf = utraj.eval(tf+hf/2);
cf = utraj.eval(tf+hf/2);

% Reconstruct SL foot Traj
xfootf = zeros(numel(pfSL), numel(tf));
for i = 1:numel(tf)
    qfk = xf(1:SL.getNumPositions(), i);
    vfk = xf(SL.getNumPositions()+(1:SL.getNumVelocities()), i);      
    kinsoli = SL.doKinematics(qfk, vfk);
    xfootf(:,i) = SL.forwardKin(kinsoli, SL.findLinkId('FLL4'), pfSL);
end

% Plot
figure(1); clf;
subplot(2,2,1); hold on; ylabel('Lift Angle')
plot(tt, xx(7,:));
plot(tf, xf(3,:), '--');
subplot(2,2,2); hold on; ylabel('Swing Angle')
plot(tt, xx(8,:));
plot(tf, xf(7,:), '--');
subplot(2,2,3); hold on; ylabel('Lift Velocity')
plot(tt+h/2, xx(nq+7,:));
plot(tf, vf(SL.getNumPositions()+1,:), '--');
subplot(2,2,4); hold on; ylabel('SwingVelocity')
plot(tt+h/2, xx(nq+8,:));
plot(tf, vf(SL.getNumPositions()+7,:), '--');

figure(2); clf;
subplot(2,2,1); hold on; ylabel('Lift Force')
yyaxis left; plot(tt+h/2, uu(1,:));
yyaxis right; plot(tf+h/2, uf(1,:), '--');
subplot(2,2,2); hold on; ylabel('Swing Force')
yyaxis left; plot(tt+h/2, uu(2,:));
yyaxis right; plot(tf+h/2, uf(2,:), '--');
subplot(2,2,3); hold on; ylabel('Normal Force')
plot(tt+h/2, cc(1,:));
plot(tf+h/2, cf(SL.nu+SL.nl+(1:SL.nc),:), '--');
subplot(2,2,4); hold on; ylabel('Tangential Basis')
plot(tt+h/2, bb(1:4,:), 'b');
plot(tf+h/2, cf(SL.nu+SL.nl+SL.nc+(1:SL.nc*SL.nd),:), 'r--');

figure(3); clf;
subplot(3,1,1); hold on; ylabel('Foot X')
plot(tt, xfoot_body(1,:));
plot(tf, xfootf(1,:), '--');

subplot(3,1,2); hold on; ylabel('Foot Y')
plot(tt, xfoot_body(2,:));
plot(tf, xfootf(2,:), '--');

subplot(3,1,3); hold on;  ylabel('Foot Z')
plot(tt, xfoot_body(3,:));
plot(tf, xfootf(3,:), '--');

