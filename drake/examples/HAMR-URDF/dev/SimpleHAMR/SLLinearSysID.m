% clear; clc; close all;
%
% urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'dev', 'SimpleHAMR', 'urdf', ...
%     'FL_scaled.urdf');
%
% % options
% options.ignore_self_collisions = true;
% options.collision_meshes = false;
% options.use_bullet = false;
% options.floating = false;
%
% SL = RigidBodyManipulator(urdf, options);
% SL = SL.setGravity([0; 0; -9.81e-3]);
% SL = compile(SL);
%
%
% nq = SL.getNumPositions();
% nv = SL.getNumVelocities();
% SL = SL.setJointLimits(-Inf(nq, 1), Inf(nq, 1));
%
% dt = 1;
% Fs = 1/dt;
%
% SL_TS = TimeSteppingRigidBodyManipulator(SL, dt/10, options);
%
% T = 6000;
% t = 0:dt:T;
% N = numel(t);
%
% f1 = 150e-3;
% band = [0, 2*f1*dt];
% uscale = 0.03;
% u = uscale*idinput([N, SL_TS.getNumInputs()], 'rgs', band)';
%
% figure(1); clf;
% subplot(2,1,1); hold on;
% plot(t, u(1,:))
% subplot(2,1,2); hold on;
% plot(t, u(2,:))
%
% freq = [-(Fs/2):(1/N/dt):(Fs/2)];
%
% f0 = 0;
% ind0 = find(freq > f0, 1, 'first');
% ind1 = find(freq > f1, 1, 'first');
% ufft = abs(fftshift((fft(u, [], 2)/size(u,2))));
%
% figure(2); clf;
% subplot(2,1,1);  hold on;
% plot(freq(ind0:ind1)*1e3, 2*ufft(1,ind0:ind1))
% subplot(2,1,2); hold on;
% plot(freq(ind0:ind1)*1e3, 2*ufft(2,ind0:ind1))
%
% %% Simulate
%
% utraj = PPTrajectory(zoh(t, u));
% utraj = utraj.setOutputFrame(SL_TS.getInputFrame);
% SL_TS_OL = cascade(utraj, SL_TS);
%
% q0 = zeros(nq,1);
% x0 = [q0; 0*q0];
% xtraj = simulate(SL_TS_OL, t, x0);
%
% tsim = xtraj.getBreaks();
% x = xtraj.eval(tsim);
%
% %%  Plot Response
% sfb_inputs = [3; 7];
% y = x(sfb_inputs, :);
% yi = interp1(tsim', y', t')';
%
% figure(1); clf;
% subplot(2,1,1); hold on;
% plot(t, u(1,:))
% plot(tsim, y(1,:));
% plot(t, rad2deg(yi(1,:)));
%
% subplot(2,1,2); hold on;
% plot(t, u(2,:))
% plot(tsim, y(2,:));
plot(t, rad2deg(yi(2,:)));

%% Compute torques about legs
xi = interp1(tsim', x', t')';
xl = zeros(2, N);
fl = zeros(2, N);

% dof = (1:SL_TS.getNumPositions())'; 
actuated_dof = SL_TS.getActuatedJoints();  
% Build stiffness
for i = 1:N
    q = xi(1:nq, i);
    qd = xi(nq+(1:nv), i);
    kinsol = SL_TS.doKinematics(q, qd, struct('compute_gradients', true));
    [K, dK] = SL_TS.positionConstraints(q);
    dKUnique = dK([1:2; 8:9; 13, 15], :); 
    Jc = -dKUnique(:,[2:4, 6:8])\dKUnique(:, actuated_dof); 
    Jc_phi = Jc([2, 5], :); 
    fl(:,i) = Jc_phi\u(:,i); 
    
%     [xli, Jli] = SL_TS.forwardKin(kinsol, SL_TS.findLinkId('FLL4'), ...
%         [0; 7.58; -11.35] , struct('rotation_type', 1));
%     xl(:,i) = xli([4, 6]);
    
%     [H, ~, B] = SL_TS.manipulatorDynamics(q, qd);
    %     fl(:,i) = J
%     qdd = SL_TS.getManipulator().dynamics(t, [q; qd], u(:,i)); % generalized acceleration
%     fl(:,i) = pinv(Jli([4,6],:)')*H*qdd(nq+(1:nv));
%     ti2 = (Jli([4,6],:)')\(H*qdd(nq+(1:nv)))
    
end

% figure(2); clf;
% subplot(2,3,1);plot(t, xl(1,:))
% subplot(2,3,2);plot(t, xl(2,:))
% subplot(2,3,3);plot(t, xl(3,:))
% 
% subplot(2,3,4); hold on;
% plot(t, rad2deg(xl(4,:))); plot(t, rad2deg(yi(2,:)));
% subplot(2,3,5);plot(t, rad2deg(xl(5,:)))
% subplot(2,3,6); hold on;
% plot(t, rad2deg(xl(6,:))); plot(t, rad2deg(yi(1,:)))


%% Plot Response Freq Domain

Fs = 1/dt;
freq = [-(Fs/2):(1/N*dt):(Fs/2)];
yfft = abs(fftshift((fft(fl, [], 2)/size(y,2))));

figure(3); clf;
subplot(2,1,1);  hold on; title('Swing')
yyaxis left; plot(freq(ind0:ind1)*1e3, 2*ufft(1,ind0:ind1))
yyaxis right; plot(freq(ind0:ind1)*1e3, 2*yfft(1,ind0:ind1))

subplot(2,1,2); hold on; title('Lift')
yyaxis left; plot(freq(ind0:ind1)*1e3, 2*ufft(2,ind0:ind1))
yyaxis right; plot(freq(ind0:ind1)*1e3, 2*yfft(2,ind0:ind1))


%% Fit linear system

inds = 1e3;
data = iddata(yi(:,inds:end)', u(:, inds:end)', dt);
ss = n4sid(data, [], 'DisturbanceModel', 'none');
[wn, z] = damp(ss);
figure(4);
compare(data, ss);


%%
v = SL_TS.constructVisualizer;
xtraj_scaled = DTTrajectory(tsim*1e-3, x);
xtraj_scaled = xtraj_scaled.setOutputFrame(xtraj.getOutputFrame());
options.slider = true;
%     xtraj.tt = xtraj.tt/1000;
v.playback(xtraj_scaled, options);

