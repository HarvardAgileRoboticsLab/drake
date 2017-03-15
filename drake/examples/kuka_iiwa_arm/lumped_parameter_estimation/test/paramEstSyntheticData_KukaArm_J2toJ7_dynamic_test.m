function paramEstSyntheticData_KukaArm_J2toJ7_dynamic_test
clear all

tmp = addpathTemporary(fullfile(pwd,'..'));

%% read arm data
path = '~/kuka-dev-estimation/drake/drake/examples/kuka_iiwa_arm/experiment_data/dynamic_test/J2toJ4_chirp';
read_joint_status_file;

joint_begin_index = 2;
joint_end_index = 7;

beginning_index = 1;
ending_index = length(joint_position_measured_raw);

index = [];
i = 1;
while i < ending_index
    if (floor(i/30) ~= floor((i+1)/30))
        index = [index, i];
    end
    i = i + 1;
end

joint_velocity_measured = joint_velocity_measured_raw(joint_begin_index:joint_end_index,index);
joint_position_measured = joint_position_measured_raw(joint_begin_index:joint_end_index,index);
joint_acceleration_measured = joint_acceleration_measured_raw(joint_begin_index:joint_end_index,index);
joint_torque_measured = joint_torque_measured_raw(joint_begin_index:joint_end_index,index);
joint_torque_measured = -joint_torque_measured;

span = 0.005;
for i=1:joint_end_index-joint_begin_index+1
    joint_acceleration_measured(i,:) = smooth(joint_acceleration_measured(i,:),span,'loess');
    joint_velocity_measured(i,:) = smooth(joint_velocity_measured(i,:),span,'loess');
    joint_torque_measured(i,:) = smooth(joint_torque_measured(i,:),span,'loess');
end

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

% Joint set to be identified
% 'J5J6J7' or 'J2J3J4'
parameterEstimationOptions.joint_set = 'J2J3J4';

% Option to print from estimator
% 'noprint'     = Do not print output from parameterEstimation.m, but will 
%                 still print output of paramEstSyntheticData.m script
% 'printEst'    = only print estimated from parameterEstimation.m
% 'printAll'	= print estimated and original from parameterEstimation.m
parameterEstimationOptions.print_result = 'noprint';

r = KukaArmPlant_J2toJ7_dynamic_test;
rCAD = KukaArmPlant_J2toJ7_dynamic_test;

outputFrameNames = r.getOutputFrame.getCoordinateNames();

arm_state = [joint_position_measured',joint_velocity_measured',joint_acceleration_measured'];
arm_input = joint_torque_measured';

outputFrameNames = [outputFrameNames;'theta2doubledot';'theta3doubledot';'theta4doubledot';'theta5doubledot';'theta6doubledot';'theta7doubledot'];
Ts = .01;
data = iddata(arm_state,arm_input,Ts,'InputName',r.getInputFrame.getCoordinateNames(),'OutputName',outputFrameNames);
[estimated_parameters, estimated_delay] = parameterEstimation_Kuka_dynamic_test(r,data,parameterEstimationOptions);

%% Print out results
coords = getCoordinateNames(r.getParamFrame);
p_true = double(rCAD.getParams);
p_init = double(r.getParams);
fprintf('\nParameter estimation results:\n\n');
fprintf('  Param  \tKuka CAD Model    \tEstimated\n');
fprintf('  -----  \t--------\t\t---------\n');
for i=1:length(coords)
  fprintf('%7s  \t%8.6f\t\t%8.6f\n',coords{i},p_true(i),estimated_parameters(i));
end