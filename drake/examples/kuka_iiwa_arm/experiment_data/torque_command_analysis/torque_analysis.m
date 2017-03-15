%=================================================
%  Date:        February 15th, 2017
%  File name:   torque_analysis.m
%  Author:      Ye Zhao (HarvardAgileRoboticsLab)
%=================================================

clear all
clc 
close all
fclose all

%% set the joint number to be identified
joint_index = 5;
%%

path = '~/kuka-dev-estimation/drake/drake/examples/kuka_iiwa_arm/experiment_data/torque_command_analysis';

read_status_file;

% crop the beginning and ending part of data sequence
beginning_index = 1;
ending_index = length(feedforward_inverse_dynamics_command);

% figure font 
bigTextSize = 20;
mediumTextSize = 18;
textSize = 16;
smallTextSize = 14;
legendMargin = 0.3;

figure(1)
title('Joints 1-4 torque components','fontsize',bigTextSize)
for i=1:4
subplot(4,1,i)
plot(feedforward_inverse_dynamics_command(i,:),'b-');
hold on;
plot(PD_impedance_ctrl_command(i,:),'r-');
hold on;
plot(total_torque_command(i,:),'y-');
xlabel('torque [Nm]','fontsize',mediumTextSize)
ylabel('time index','fontsize',mediumTextSize)
grid on
end
% xlim([-1.2, 1.2])
% ylim([-0, 0.2])
print -depsc Joint_1-4_torque_components

figure(2)
title('Joints 5-7 torque components','fontsize',bigTextSize)
for i=5:7
subplot(3,1,i-4)
plot(feedforward_inverse_dynamics_command(i,:),'b-');
hold on;
plot(PD_impedance_ctrl_command(i,:),'r-');
hold on;
plot(total_torque_command(i,:),'y-');
xlabel('time index','fontsize',mediumTextSize)
ylabel('joint torque [Nm]','fontsize',mediumTextSize)
grid on
end
% xlim([-1.2, 1.2])
% ylim([-0, 0.2])
print -depsc Joint_5-7_torque_components

