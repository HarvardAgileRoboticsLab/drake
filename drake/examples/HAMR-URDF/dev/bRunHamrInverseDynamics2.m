clear; close all; clc;

%% Build Robot

hamr_urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'dev', 'SimpleHAMR', ...
    'urdf', 'HAMRSimple_scaled.urdf');

% options
options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = true;       % CHANGE
options.collision = true;

hS = HAMRSimpleRBM(hamr_urdf,options);
x0 = hS.getInitialState();

nq = hS.getNumPositions();
nv = hS.getNumVelocities();
nu = hS.getNumInputs();
nc = hS.getNumContactPairs();
nd = 4; % pyramidal friction cone approx

%% Load Trajectory

fname = 'TrajOpt_MovingBody_SimpleSprings7';
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

%% Form foot positions

leg_links = {'FL2', 'RL2', 'FR2', 'RR2'};     % leg links

pfhamr = [0 0 -14.988382167532292]'; % position of foot in local frame

xfoot = zeros(numel(pfhamr)*nc, numel(tt));
for i = 1:numel(tt)
    kinsol = hS.doKinematics(xx(1:nq, i), zeros(nv,1));
    for j = 1:nc
        xfoot(j*numel(pfhamr)-2:j*numel(pfhamr),i) = ...
            hS.forwardKin(kinsol, hS.findLinkId(leg_links{j}), pfhamr);
    end
end

%% Solve

urdf = 'HAMR_scaledV2';
feet = {'FLL4', 'RLL4', 'FRL4', 'RRL4'};
feet_pos = [0, 7.58, -11.35;
    0, 7.58, -11.35;
    0, -7.58, -11.35;
    0, -7.58, -11.35]';

td.xfoot_traj = PPTrajectory(foh(tt, xfoot));
tic;
[hamr,xtraj,utraj,z,F,info] = HamrInvDynDirtran(tt, urdf, feet, feet_pos, td, options);
toc

%% Plotting and playback

% Evaluate final trajectory
tf = xtraj.getBreaks();
hf = mean(diff(tf));
xf = xtraj.eval(tf);
vf = xtraj.eval(tf+hf/2);
uf = utraj.eval(tf+hf/2);
cf = utraj.eval(tf+hf/2);

% Swing and lift angles
SFB_INPUTS = [3, 7; 11, 15; 19, 23; 27, 31];
trans = {'Front Left', 'Rear Left', 'Front Right', 'Rear Right'};

% Floating base
if options.floating
    nfb = 6;
else
    nfb = 0;
end

% Reconstruct SL foot Traj
xfootf = zeros(numel(feet_pos(:)), numel(tf));
for i = 1:numel(tf)
    qfk = xf(1:hamr.getNumPositions(), i);
    vfk = 0*qfk; 
    kinsol = hamr.doKinematics(qfk, vfk);
    for j = 1:nc
        xfootf(j*numel(pfhamr)-2:j*numel(pfhamr),i) = ...
            hamr.forwardKin(kinsol, hamr.findLinkId(feet{j}), feet_pos(:,j));
    end
end

% Plot Swing and lift angles
figure(1); clf;
for j = 1:nc
    subplot(2,nc,j); hold on;
    title(trans{j}); ylabel('Swing Angle')
    plot(tt, rad2deg(xx(nfb+2*j-1,:)));
    plot(tf, rad2deg(xf(nfb+SFB_INPUTS(j,1),:)), '--');
    
    subplot(2,nc,nc+j); hold on;
    ylabel('Lift Angle'); xlabel('Time(ms)')
    plot(tt, rad2deg(xx(nfb+2*j,:) - x0(nfb+2*j)));
    plot(tf, rad2deg(xf(nfb+SFB_INPUTS(j,2),:)), '--');
end

% Plot input forces
figure(2); clf;
for j = 1:nc
    subplot(2,nc,j); hold on;
    title(trans{j}); ylabel('Swing Force')
    plot(tt+h/2, uu(2*j-1,:));
    plot(tf+h/2, uf(2*j-1,:), '--');
    
    subplot(2,nc,nc+j); hold on;
    ylabel('Lift Force'); xlabel('Time(ms)')
    plot(tt+h/2, uu(2*j,:));
    plot(tf+h/2, uf(2*j-1,:), '--');
end


figure(3); clf; hold on;
for j = 1:nc
subplot(4,1,j); hold on; 
ylabel('Normal Force'); xlabel('Time(ms)')
title(feet{j}); 
plot(tt+h/2, cc(j,:));
plot(tf+h/2, cf(hamr.nu+hamr.nl+j,:), '--');
end
% Plot foot positions
figure(4); clf;
for j = 1:nc
    subplot(3,nc,j); hold on; 
    title(feet{j});ylabel('Foot X')
    plot(tt, xfoot((j-1)*numel(pfhamr)+1,:));
    plot(tf, xfootf((j-1)*numel(pfhamr)+1,:), '--');
    
    subplot(3,nc,j+nc); hold on; 
    ylabel('Foot Y')
    plot(tt, xfoot((j-1)*numel(pfhamr)+2,:));
    plot(tf, xfootf((j-1)*numel(pfhamr)+2,:), '--');
    
    subplot(3,nc,j+2*nc); hold on;  
    ylabel('Foot Z'); xlabel('Time(ms)')
    plot(tt, xfoot((j-1)*numel(pfhamr)+3,:));
    plot(tf, xfootf((j-1)*numel(pfhamr)+3,:), '--');
end

if options.floating
    figure(5); clf;
    title_str = {'Com-X(mm)', 'Com-Y(mm)', 'Com-Z(mm)', 'Roll(deg)', 'Pitch(deg)', 'Yaw(deg)'};
    for i = 1:6
        subplot(2,3,i); hold on; ylabel(title_str(i), 'FontSize', 18)
        xlabel('Time(ms)', 'FontSize', 18)
        if i > 3
            plot(tt, rad2deg(xx(i,:)));             
            plot(tf, rad2deg(xf(i,:)), '--');            
        else
            plot(tt, xx(i,:)); 
            plot(tf, xf(i,:), '--');                                     
        end
    end    
end

%% Playback and save
xtraj_scaled = DTTrajectory(tf*1e-3, xf);
xtraj_scaled = xtraj_scaled.setOutputFrame(xtraj.getOutputFrame());
options.slider = true;
v = hamr.constructVisualizer(); 
v.playback(xtraj_scaled, options);
save([fname, '_DynfullRobot'], 'xtraj', 'utraj')