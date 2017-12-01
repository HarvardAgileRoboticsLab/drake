clear; clc; close all; 
fname = 'TrajOpt_MovingBody_SimpleSprings7'; 

%% Build Robot 

% file
urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'dev', 'SimpleHAMR', 'urdf', 'HAMRSimple_scaled.urdf');

% options
options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = true;
options.collision = true;

hamr = HAMRSimpleRBM(urdf,options);
x0 = hamr.getInitialState();
v = hamr.constructVisualizer();

%%  Plot

load(fname); 

nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nx = nq+nv;
nu = hamr.getNumInputs();

tt = xtraj.getBreaks();
h = mean(diff(tt));
xx = xtraj.eval(tt);
vv = xtraj.eval(tt); 
uu = utraj.eval(tt);
ss = straj.eval(tt); 

act_dof = hamr.getActuatedJoints();
nq = hamr.getNumPositions();
title_str = {'Front Left Swing', 'Front Left Lift', ...
    'Rear Left Swing', 'Rear Left Lift', ...
    'Front Right Swing', 'Front Right Lift', ...
    'Rear Rear Swing', 'Rear Rear Lift'};

figure(1); clf; hold on;
for i = 1:numel(act_dof)
    subplot(4,2,i); hold on; title(title_str(i))
    yyaxis left; hold on; plot(tt, uu(i,:)); ylabel('Force(N)')
    yyaxis right; hold on; plot(tt, rad2deg(xx(act_dof(i), :)-x0(act_dof(i)))); 
    ylabel('Deflection(deg)')
%     legend('Force', 'Deflection')
end

figure(2); clf;
title_str = {'com x', 'com y', 'com z', 'roll', 'pitch', 'yaw'};
for i = 1:6
    subplot(2,3,i); hold on; title(title_str(i))
    yyaxis left; hold on; plot(tt, xx(i,:)); %ylabel('Position')
    yyaxis right; hold on; plot(tt, xx(i+nq, :)); %ylabel('Velocity')
    
    if i > 3
        yyaxis left; hold on; plot(tt, rad2deg(xx(i,:))); %ylabel('Position')
        yyaxis right; hold on; plot(tt, rad2deg(xx(i+nq, :))); %ylabel('Velocity')
    end
%     legend('Position', 'Velocity')
end

angle_inds = [hamr.body(~[hamr.body.floating]).position_num]';
angle_inds = angle_inds(angle_inds>0);

phi = zeros(4, numel(tt));
cc = ctraj.eval(tt + h/2);
for i = 2:numel(tt)
    q = xx(1:nq, i);
%     q2 = xx(1:nq, i+1);
    qd = vv(nq+1:nx, i);
    kinsol = hamr.doKinematics(q, qd);
    phi(:,i-1) = hamr.contactConstraints(kinsol, false);
end


legs = {'FL2', 'RL2', 'FR2', 'RR2'};

figure(3); clf; hold on; 
for i = 1:size(phi, 1)    
    subplot(2,2,i); hold on; title(legs{i}); 
    yyaxis left; plot(tt+h/2, phi(i,:), '*-'); ylabel('Z-Distance(mm)'); %ylim([0, 5])
    yyaxis right; plot(tt+h/2, cc(i,:), '*-'); ylabel('Force (N)'); %ylim([0, 20e-3])
end

figure(4); clf; hold on; 
subplot(2,1,1); hold on; 
for i = 1:numel(tt)
    nc(i) = phi(:,i)'*cc(:,i); 
end
yyaxis left; plot(tt, nc);
yyaxis right; plot(tt, ss, '--')
legend('Normal Constraint', 'S')
subplot(2,1,2); hold on; 
plot(tt, nc'-ss');
plot(tt, 1e-5*ones(numel(tt), 1), 'k'); 


lp_b = [0, 0, -14.97;
    0, 0, -14.97;
    0, 0, -14.97;
    0, 0, -14.97];

lp_g = zeros([numel(tt), size(lp_b')]);

legs = {'FL2', 'RL2', 'FR2', 'RR2'};

for j = 1:numel(tt)
    q = xx(1:nq, j);
    qd = vv(nq+1:nx, j);
    kinsol = hamr.doKinematics(q, qd);
    for i = 1:size(lp_b,1)
        lp_g(j,:,i) = hamr.forwardKin(kinsol, hamr.findLinkId(legs{i}), lp_b(i,:)');
    end
end

figure(5); clf; hold on;
% figure(j); clf;
subplot(3,1,1); hold on; ylabel('Foot X')
for i = 1:size(lp_b,1)
    plot(tt, lp_g(:,1,i))
end
subplot(3,1,2); hold on; ylabel('Foot Y')
for i = 1:size(lp_b,1)
    plot(tt, lp_g(:,2,i))
end

subplot(3,1,3); hold on;  ylabel('Foot Z')
for i = 1:size(lp_b,1)
    plot(tt, lp_g(:,3,i))
end
legend(legs)

%% Animate 
qq = xtraj.eval(tt); qq = qq(1:hamr.getNumPositions(), :); 
% qq1 = qq; 
% qq1(7:2:13, :) = qq(8:2:14, :); 
% qq1(8:2:14, :) = qq(7:2:13, :); 
% qq = qq1; 
qtraj_scaled = PPTrajectory(foh(tt*1e-3, qq));
qtraj_scaled = qtraj_scaled.setOutputFrame(v.getInputFrame());
v.playback(qtraj_scaled, struct('slider', true));
