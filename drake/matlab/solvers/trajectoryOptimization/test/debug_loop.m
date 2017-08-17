clear; clc; close all;

% Plant
options.floating = false;
options.use_bullet = false;
file = fullfile(getDrakePath,'examples', 'SimpleFourBar', 'FourBar_JointLimits4DOF.urdf');
plant = VariationalFourBarPlant(file, options);
plant = compile(plant);

v = plant.constructVisualizer(); 

% Time% Steps
N=31;

% Initial condition
q0 = pi/2*ones(4,1); %[9*pi/10; pi/10; 9*pi/10; pi/10]; %pi/4*ones(3,1);
valuecheck(positionConstraints(plant, q0), zeros(6,1), 1e-4);
x0 = [q0;0*q0];
% v.inspector(x0); 

% Final time
tf = 4;

% simtraj = plant.simulate([0, tf], x0); 
% v.playback(simtraj, struct('slider', true))

%% Simulates

[q, kl] = variationalFourBarForwardSim(plant, N, x0, tf);
[qopt,utraj,ctraj,btraj,psitraj,etatraj,jltraj, klopt, straj] = variationalFourBar(plant, N, x0, tf);
%%

qsim = q.eval(q.getBreaks());
qqopt = qopt.eval(qopt.getBreaks());

figure(3); clf;
for i=1:size(q, 1)
    subplot(2,size(q,1)/2,i); hold on;
    plot(rad2deg(qsim(i,:)), 'b');
    if i == 1
        plot(rad2deg(10*pi/11)*ones(size(qsim(i,:))), 'k'); 
    end
    plot(rad2deg(qqopt(i,:)), 'r--');
    hold off;
    legend('MidpointRule', 'TrajOpt')
end

% good_ind= [1; 3];
kllopt = circshift(klopt.eval(klopt.getBreaks()),1,2);
% kllopt = kllopt(good_ind, :);
figure(4); clf;
for i=1:size(kllopt, 1)
    subplot(size(kllopt,1),1,i); hold on;
    if i == 1
        plot(kl(i,1:end), 'b');
    end
    if i == 2
        plot(kl(2,1:end), 'b');
    end
    if i == 3
        plot(kl(3,1:end), 'b');
    end
    plot(kllopt(i,1:end), 'r--');
    hold off;
    legend('MidpointRule', 'TrajOpt')
end
tilefigs; 

v.playback(qopt, struct('slider', true)); 
