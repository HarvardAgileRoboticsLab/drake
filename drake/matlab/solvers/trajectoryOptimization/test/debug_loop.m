clear; clc; close all;

% Plant
options.floating = false;
options.use_bullet = false;
file = fullfile(getDrakePath,'examples', 'SimpleFourBar', 'FourBar_JointLimits.urdf');
plant = RigidBodyManipulator(file,options);

% Time Steps
N=21;

% Initial condition
q0 = [9*pi/10; pi/10; 9*pi/10]; %pi/4*ones(3,1);
valuecheck(positionConstraints(plant, q0), zeros(6,1), 1e-4);
x0 = [q0;0*q0];

% Final time
tf = 2;

%% Simulates

[q, kl] = variationalFourBarForwardSim(plant, N, x0, tf);
[qopt,utraj,ctraj,btraj,psitraj,etatraj,jltraj, klopt] = variationalFourBar(plant, N, x0, tf);
%%

qsim = q.eval(q.getBreaks());
qqopt = qopt.eval(qopt.getBreaks());

figure(2); clf;
for i=1:size(q, 1)
    subplot(2,size(q,1)/2,i); hold on;
    plot(qsim(i,:), 'b');
    plot(qqopt(i,:), 'r--');
    hold off;
    legend('MidpointRule', 'TrajOpt')
end

% good_ind= [1; 3];
kllopt = circshift(klopt.eval(klopt.getBreaks()),1,2);
% kllopt = kllopt(good_ind, :);
figure(3); clf;
for i=1:size(kllopt, 1)
    subplot(2,size(kllopt,1)/2,i); hold on;
    if i == 1
        plot(kl(i,1:end), 'b');
    end
    if i == 2
        plot(kl(2,1:end), 'b');
    end
    plot(kllopt(i,1:end), 'r--');
    hold off;
    legend('MidpointRule', 'TrajOpt')
end
