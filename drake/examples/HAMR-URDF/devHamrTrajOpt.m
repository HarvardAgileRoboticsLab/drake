clear; clc; close all;

[hamr,xtraj,utraj,ctraj,btraj,...
    psitraj,etatraj,jltraj, kltraj, straj] = runVariationalTrajOpt();

save(['TrajOpt_', date], 'xtraj', 'utraj', 'ctraj', 'btraj', 'psitraj', 'etatraj', ...
    'jltraj', 'kltraj', 'straj'); 

%% Playback

v = hamr.constructVisualizer();
xtraj_scaled = PPTrajectory(foh(xtraj.getBreaks()*1e-3, xtraj.eval(xtraj.getBreaks())));
xtraj_scaled = xtraj_scaled.setOutputFrame(xtraj.getOutputFrame());
v.playback(xtraj_scaled, struct('slider', true));


%% Plotting
% load('NomTraj.mat')
% ttnom = tt; 
% uunom = yy(101:108,:);

% tt = xtraj.getBreaks();
% yy = xtraj.eval(tt);
% uu = utraj.eval(tt);

act_dof = hamr.getActuatedJoints();
nq = hamr.getNumPositions();
title_str = {'Front Left Swing', 'Front Left Lift', ...
    'Rear Left Swing', 'Rear Left Lift', ...
    'Front Right Swing', 'Front Right Lift', ...
    'Rear Rear Swing', 'Rear Rear Lift'};

figure(1); clf; hold on;
for i = 1:numel(act_dof)
    subplot(4,2,i); hold on; title(title_str(i))
    plot(tt, uu(i, :)*1e3); hold on;
    plot(ttnom, uunom(i,:)*1e3);   
    %yyaxis right; plot(tt, yy(act_dof(i), :)*1e3)
    %legend('Force(mN)', 'Deflection(\mum)')
end

figure(2); clf;
title_str = {'com x', 'com y', 'com z', 'roll', 'pitch', 'yaw'};
for i = 1:6
    subplot(3,2,i); hold on; title(title_str(i))
    yyaxis left; hold on; plot(tt, yy(i,:))
    yyaxis right; hold on; plot(tt, yy(i+nq, :))
    legend('Position', 'Velocity')
end

% phi = zeros(4, numel(tt));
% for i = 1:numel(tt)
%     q = yy(1:ndof/2, i);
%     qd = yy(ndof/2+1:ndof, i);
%     kinsol = hamr.doKinematics(q, qd);
%     phi(:,i) = hamr.contactConstraints(kinsol, false, contact_opt);
% end
% 
% figure(3); clf; hold on;
% plot(tt, phi);