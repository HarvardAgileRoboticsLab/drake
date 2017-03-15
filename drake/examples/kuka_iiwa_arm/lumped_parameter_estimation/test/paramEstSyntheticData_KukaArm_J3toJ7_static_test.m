function paramEstSyntheticData_KukaArm_J3toJ7
clear all

tmp = addpathTemporary(fullfile(pwd,'..'));

%% read arm data
path = '~/kuka-dev-estimation/drake/drake/examples/kuka_iiwa_arm/experiment_data/static_test/J5toJ7_set2';
read_joint_status_file;

% % crop the beginning and ending part of data sequence
% beginning_index_set = 560;%[560,2530,7270];%[3600, 6860];%;%2800;%1000;%
% ending_index_set = 10000;%length(joint_position_measured_raw) - 400;%[1960,4360,8980];%[6200, 11000];%;%;%1800;%%1800;%1800;%3600;%2000;%
% 
% % bad index
% % 4900,
% % 6740,
% 
% joint_velocity_measured = [];
% joint_position_measured = [];
% joint_acceleration_measured = [];
% joint_torque_measured = [];
% 
% for i=1:length(beginning_index_set)
%     joint_velocity_measured = [joint_velocity_measured, joint_velocity_measured_raw(3:7,beginning_index_set(i):ending_index_set(i))];
%     joint_position_measured = [joint_position_measured, joint_position_measured_raw(3:7,beginning_index_set(i):ending_index_set(i))];
%     joint_acceleration_measured = [joint_acceleration_measured, joint_acceleration_measured_raw(3:7,beginning_index_set(i):ending_index_set(i))];
%     joint_torque_measured = [joint_torque_measured, joint_torque_measured_raw(3:7,beginning_index_set(i):ending_index_set(i))];
% end
% 
% joint_torque_measured = -joint_torque_measured;

joint_begin_index = 3;
joint_end_index = 7;

beginning_index = 1;
ending_index = length(joint_position_measured_raw);

joint_velocity_measured = joint_velocity_measured_raw(joint_begin_index:joint_end_index,beginning_index:ending_index);
joint_position_measured = joint_position_measured_raw(joint_begin_index:joint_end_index,beginning_index:ending_index);
joint_acceleration_measured = joint_acceleration_measured_raw(joint_begin_index:joint_end_index,beginning_index:ending_index);
joint_torque_measured = joint_torque_measured_raw(joint_begin_index:joint_end_index,beginning_index:ending_index);
joint_torque_measured = -joint_torque_measured;

% q5 = joint_position_measured(1,:);
% q6 = joint_position_measured(2,:);
% q7 = joint_position_measured(3,:);
% 
% dq5 = joint_velocity_measured(1,:);
% dq6 = joint_velocity_measured(2,:);
% dq7 = joint_velocity_measured(3,:);
% 
% q5ddot = joint_acceleration_measured(1,:);
% q6ddot = joint_acceleration_measured(2,:);
% q7ddot = joint_acceleration_measured(3,:);
% 
% 
%     l5x = 0, l5y = 0.1845, l5z = 0;
%     l6x = 0, l6y = 0, l6z = 0.2155;
%     l7x = 0,l7y = 0.081, l7z = 0;
%     
%     %m1 = 5.76; m2 = 6.35; m3 = 3.5; m4 = 3.5; 
%     m5 = 3.5; m7 = 1.2;
%     g = 9.81;
%     
%     c5x = 0.0001, c5y = 0.021, c5z = 0.076;
%     c6x = 0, c6y = 0.0006, c6z = 0.0004;
%     c7x = 0, c7y = 0, c7z = 0.02;
%     
% %     I1xx= 0.033, I1xy= 0, I1xz= 0, I1yy= 0.0333, I1yz= 0.004887, I1zz= 0.0123;
% %     I2xx= 0.0305, I2xy= 0, I2xz= 0, I2yy= 0.0304, I2yz= 0.004887, I2zz= 0.011;
% %     I3xx= 0.025, I3xy= 0, I3xz= 0, I3yy= 0.0238, I3yz= 0.00487, I3zz= 0.0076;
% %     I4xx= 0.017, I4xy= 0, I4xz= 0, I4yy= 0.0164, I4yz= 0.00284, I4zz= 0.006;
%     I5xx= 0.01, I5xy= 0, I5xz= 0, I5yy= 0.0087, I5yz= 0.00309, I5zz= 0.00449;
%     I6xx= 0.0049, I6xy= 0, I6xz= 0, I6yy= 0.0047, I6yz= 0.000246, I6zz= 0.0036;
%     I7xx= 0.0002, I7xy= 0, I7xz= 0, I7yy= 0.0002, I7yz= 0, I7zz= 0.0003;
%     
% w = (- c6y*c6z*cos(q6).*q5ddot) + c6y^2 .*q6ddot - (c6y^2*dq5.^2.*sin(2*q6))/2 - g*sin(q6)*c6y;
% 
% H1 = (2268949521066275.*cos(q6))./9223372036854775808 - I7xx*cos(q7).*sin(q6).*sin(q7) + I7yy*cos(q7).*sin(q6).*sin(q7).*q5ddot;
% H2 = (m7*c7z^2 + (81*m7*c7z)/500 + I7xx/2 + I7yy/2 + I6zz + (6561*m7)/1000000 - (I7xx*cos(2*q7))/2 + (I7yy*cos(2*q7))/2).*q6ddot;
% 
% C2 = dq5.*((I7xx*dq6.*sin(2*q5))/2 - (I7yy*dq6.*sin(2*q5))/2 + I7yy*dq7.*cos(q5).^2.*sin(q6) + I7xx*dq7.*sin(q5).^2.*sin(q6)) - (dq5.^2*m7.*sin(2*q6)*(1000*c7z + 81)^2)/2000000 + (dq5.*dq6.*sin(2*q5)*(I6xx - I6yy))/2 - dq6.*dq7.*cos(q5).*cos(q6).*sin(q5)*(I7xx - I7yy);
% 
% G2 = -(g*sin(q6)*(81*m7 + 1000*c7z*m7))/1000;
% 
% y = joint_torque_measured(2,:) - H1 - H2 - C2 - G2;
% 
% figure(1)
% plot(w,y)

%%
% Synthetic parameter estimation mode
% 'base':       Base case parameter estimation - intial params = true params
% 'paramerr':   w/ parameter error but no measurement noise
% 'measnoise':  w/ parameter error and w/ measurement noise
% 'delay':      w/ parameter error and w/ measurement noise and delay (Not complete)
mode = 'measnoise';

% Parameter Estimation model
% 'dynamic'     = use dynamic model - requires qdd
% 'energetic'   = use energetic model - doesn't require qdd
parameterEstimationOptions.model = 'dynamic';

% Parameter Estimation robot
% 'Acrobot'
% 'KukaArm'
parameterEstimationOptions.robot = 'KukaArm';

% Method by which to obtain qdd (not used if using energetic model)
% 'manipul':   Use acrobot manipulator equations to estimate true qdd
% 'derivative': Take the derivative of qd
qddmode = 'manipul';

% Parameter Estimation method
% 'nonlinprog'  = nonlinear least squares (LS) to solve problem
% 'linprog'       = linear LS on lumped params then nonlinear LS to recover
%                 original parameters
% 'simerr'      = minimize simulation error
% 'lsqnonlin'   = use MATLAB's built-in nonlinear least squares solver (debugging)
parameterEstimationOptions.method = 'nonlinprog';

% Option to print from estimator
% 'noprint'     = Do not print output from parameterEstimation.m, but will 
%                 still print output of paramEstSyntheticData.m script
% 'printEst'    = only print estimated from parameterEstimation.m
% 'printAll'	= print estimated and original from parameterEstimation.m
parameterEstimationOptions.print_result = 'noprint';

% Standard deviation of the data input error
if strcmp(parameterEstimationOptions.model,'energetic')
    % In this case, for theta1,theta2,theta1dot,theta2dot
    noisestd = sqrt([.0005, .0005, .0005, .0005, .0005, .0005, .0005, .0007, .0007, .0007, .0007, .0007, .0007, .0007]);
else
    % In this case, for theta1,theta2,theta1dot,theta2dot,theta1doubledot,theta2doubledot
    noisestd = sqrt([.0005, .0005, .0005, .0005, .0005, .0005, .0005, .0007, .0007, .0007, .0007, .0007, .0007, .0007, .0012, .0012, .0012, .0012, .0012, .0012, .0012]);
end
% Standard deviation of the parameter value percent error
paramstd = 1/5;

r = KukaArmPlant_J3toJ7;
rtrue = KukaArmPlant_J3toJ7;

% if ~strcmp(mode,'base')
%     % Perturb original parameter estimates with random percentage error
%     % normally distributed with standard dev = paramstd, and greater than -1
%     paramerr = randn(1,10)*paramstd;
%     while sum(paramerr<=-1)~=0
%         paramerr(paramerr<-1) = randn(1,sum(paramerr<-1))*paramstd;
%     end
% 
%     % rtrue.l1 = rtrue.l1 + rtrue.l1*paramerr(1); 
%     % rtrue.l2 = rtrue.l2 + rtrue.l2*paramerr(2); 
%     % rtrue.m1 = rtrue.m1 + rtrue.m1*paramerr(3); 
%     % rtrue.m2 = rtrue.m2 + rtrue.m2*paramerr(4);
%     rtrue.bv1_positive  = rtrue.bv1_positive + rtrue.bv1_positive*paramerr(5);
%     rtrue.bv2_positive  = rtrue.bv2_positive + rtrue.bv2_positive*paramerr(6);
%     rtrue.c1x = rtrue.c1x + rtrue.c1x*paramerr(7); 
%     rtrue.c2x = rtrue.c2x + rtrue.c2x*paramerr(8); 
%     rtrue.I1xx = rtrue.I1xx + rtrue.I1xx*paramerr(9);  
%     rtrue.I2xx = rtrue.I2xx + rtrue.I2xx*paramerr(10);
% end

outputFrameNames = r.getOutputFrame.getCoordinateNames();

% %% Test on swingup up data
% % [utraj,xtraj] = swingUpTrajectory(rtrue);
% % Ts = .01; breaks=getBreaks(utraj); T0 = breaks(1); Tf = breaks(end);
% % tsamples = T0:Ts:Tf;
% Ts = .01; 
% tsamples = 0:Ts:6;
% utraj = zeros(601,7);
% xtraj = zeros(601,14);
% xsamples = zeros(601,14);%eval(xtraj,tsamples)';
% usamples = zeros(601,7);%eval(utraj,tsamples)';

% % calcuate qddot
% nq = r.num_positions;
% qdd = zeros(length(tsamples),nq);
% for i=1:length(tsamples)
%     [H,C,B] = manipulatorDynamics(rtrue,xsamples(i,1:nq)',xsamples(i,nq+(1:nq))');
%     qdd(i,:) = (H\(B*usamples(i,:)' - C))';
% end

% %% Generate second derivative
% if ~strcmp(parameterEstimationOptions.model,'energetic')
%     nq = r.num_positions;
%     if strcmp(qddmode,'manipul')
%         qdd = zeros(length(tsamples),nq);
%         for i=1:length(tsamples)
%           [H,C,B] = manipulatorDynamics(rtrue,xsamples(i,1:nq)',xsamples(i,nq+(1:nq))');
%           qdd(i,:) = (H\(B*usamples(i,:)' - C))';
%         end
%     else
%         % Differentiating to get the second derivative of the state variables
%         % TODO: Try lowpass filter
%         qdd = deriv(xtraj,T0:Ts:Tf)'; qdd = qdd(:,nq+(1:nq));
%     end
%     xsamples = [xsamples qdd];
%     outputFrameNames = [outputFrameNames;'theta1doubledot';'theta2doubledot'];
%     
% %     qddTrue = qdd;
% %     qdd = deriv(xtraj,T0:Ts:Tf)'; qdd = qdd(:,nq+(1:nq));
% %     figure;plot(qddTrue(:,1),'-g'); hold on; plot(qdd(:,1),'-r')
% %     title(['q_1 True Acceleration vs. Derivative of Velocity']);
% %     xlabel('Sample Time')
% %     ylabel({'Acceleration Magnitude' '(m/s)/sample'})
% %     legend('True Acceleration','Derivative of Vel');
% end

% Only for simulated data
% %% Add gaussian noise to measurements
% if strcmp(mode,'measnoise') || strcmp(mode,'delay')
%     measurementNoise = randn(size(xsamples))*diag(noisestd);
% else
%     measurementNoise = 0;
% end
% xsamplesfinal = xsamples+measurementNoise;

arm_state = [joint_position_measured',joint_velocity_measured',joint_acceleration_measured'];
arm_input = joint_torque_measured';

% % Debugging purposes
% v = AcrobotVisualizer(r);
% vtrue = AcrobotVisualizer(rtrue);
% vtrue.playback(xtraj);
% 
% % plot(xsamples(:,1)); hold on; plot(xsamplesfinal(:,1));
% plot(xsamples); hold on; plot(xsamplesfinal);
% % Testing to see if estimate is produced by local minima
% r = rtrue;

outputFrameNames = [outputFrameNames;'theta3doubledot';'theta4doubledot';'theta5doubledot';'theta6doubledot';'theta7doubledot'];
Ts = .01;
data = iddata(arm_state,arm_input,Ts,'InputName',r.getInputFrame.getCoordinateNames(),'OutputName',outputFrameNames);
[estimated_parameters, estimated_delay] = parameterEstimation(r,data,parameterEstimationOptions);

%% Print out results
coords = getCoordinateNames(r.getParamFrame);
p_true = double(rtrue.getParams);
p_init = double(r.getParams);
fprintf('\nParameter estimation results:\n\n');
fprintf('  Param  \tKuka CAD Model    \tEstimated\n');
fprintf('  -----  \t--------\t\t---------\n');
for i=1:length(coords)
  fprintf('%7s  \t%8.6f\t\t%8.6f\n',coords{i},p_true(i),estimated_parameters(i));
end
