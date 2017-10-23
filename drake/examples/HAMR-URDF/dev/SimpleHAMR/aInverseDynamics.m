clear; close all; clc; 

%% Build Robot

hamr_urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'dev', 'SimpleHAMR', ...
    'urdf', 'HAMRSimple_scaled.urdf');

% options
options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = false;
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
uu = td.utraj.eval(tt);
cc = td.ctraj.eval(tt);
bb = td.btraj.eval(tt);

%% Form foot position and Forces

pf = [0 0 -14.97]'; % position of foot in local frame
xfoot_body = zeros(numel(pf), numel(tt));
xfoot_world = zeros(numel(pf), numel(tt));
Fext = zeros(3, numel(tt));

dkopt.compute_gradients = true;
fkopt.base_or_frame_id = hamr.findLinkId('Chassis');

for i = 1:numel(tt)
    kinsol = hamr.doKinematics(xx(1:nq, i), zeros(nv,1), dkopt);
    [xfoot_body(:,i)] = hamr.forwardKin(kinsol, hamr.findLinkId('FL2'), pf, fkopt);
end

%% Solve
xfoot_traj = PPTrajectory(foh(tt, xfoot_body));
tic;
[SL, xtraj,utraj,z,F,info] = SLInverseDynDirtran(tt, xfoot_traj, ...
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
% uf = utraj.eval(tf);
% cf = utraj.eval(tf);
% 
% % Reconstruct SL foot Traj
% xfootf = zeros(numel(pf), numel(tf));
% for i = 1:numel(tf)
%     qfk = xf(1:SL.getNumPositions(), i);
%     vfk = xf(SL.getNumPositions()+(1:SL.getNumVelocities()), i);
%     
%     kinsoli = SL.doKinematics(qfk, vfk);
%     xfootf(:,i) = SL.forwardKin(kinsoli, SL.findLinkId('FLL4'), pf);
% end

% Plot
figure(1); clf;
subplot(2,2,1); hold on; title('Lift Angle')
% plot(tt, xx(1,:));
plot(tf, xf(8,:), '--');
subplot(2,2,2); hold on; title('Swing Angle')
% plot(tt, xx(2,:));
plot(tf, xf(4,:), '--');
subplot(2,2,3); hold on; title('Lift Angular Velocity')
% plot(tt+h/2, xx(nq+1,:));
plot(tf, xf(SL.getNumPositions()+8,:), '--');
subplot(2,2,4); hold on; title('Swing Angular Velocity')
% plot(tt+h/2, xx(nq+2,:));
plot(tf, xf(SL.getNumPositions()+4,:), '--');

% figure(2); clf;
% subplot(2,2,1); hold on; title('Lift Force')
% plot(tt+h/2, uu(1,:));
% plot(tf, uf(1,:), '--');
% subplot(2,2,2); hold on; title('Swing Force')
% plot(tt+h/2, uu(2,:));
% plot(tf, uf(2,:), '--');
% subplot(2,2,3); hold on; title('Normal Force')
% plot(tt, cc(1,:));
% % plot(tf, uf(3,:), '--');
% subplot(2,2,4); hold on; title('Tangential Force Basis')
% plot(tt, bb(1:4,:), 'b');
% % plot(tf, uf(4:7,:), 'r--');

% figure(3); clf;
% subplot(3,1,1); hold on;
% plot(tt, xfoot_body(1,:));
% plot(tf, xfootf(1,:), '--');
% 
% subplot(3,1,2); hold on;
% plot(tt, xfoot_body(2,:));
% plot(tf, xfootf(2,:), '--');
% 
% subplot(3,1,3); hold on;
% plot(tt, xfoot_body(3,:));
% plot(tf, xfootf(3,:), '--');

