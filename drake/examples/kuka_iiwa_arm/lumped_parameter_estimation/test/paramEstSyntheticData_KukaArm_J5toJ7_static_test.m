function paramEstSyntheticData_KukaArm_J5toJ7_static_test
clear all

tmp = addpathTemporary(fullfile(pwd,'..'));

%% read arm data
path = '~/kuka-dev-estimation/drake/drake/examples/kuka_iiwa_arm/experiment_data/static_test/J5toJ7';
read_joint_status_file;

% extract all the joint data
joint_index = 6;
joint_begin_index = 5;
joint_end_index = 7;

beginning_index = 1;
ending_index = length(joint_position_measured_raw);

joint_velocity_measured = joint_velocity_measured_raw(joint_begin_index:joint_end_index,beginning_index:ending_index);
joint_position_measured = joint_position_measured_raw(joint_begin_index:joint_end_index,beginning_index:ending_index);
joint_acceleration_measured = joint_acceleration_measured_raw(joint_begin_index:joint_end_index,beginning_index:ending_index);
joint_torque_measured = joint_torque_measured_raw(joint_begin_index:joint_end_index,beginning_index:ending_index);
joint_torque_measured = -joint_torque_measured;

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

r = KukaArmPlant_J5toJ7_static_test;
rCAD = KukaArmCADPlant;

outputFrameNames = r.getOutputFrame.getCoordinateNames();

arm_state = [joint_position_measured',joint_velocity_measured',joint_acceleration_measured'];
arm_input = joint_torque_measured';

outputFrameNames = [outputFrameNames;'theta5doubledot';'theta6doubledot';'theta7doubledot'];
Ts = .01;
data = iddata(arm_state,arm_input,Ts,'InputName',r.getInputFrame.getCoordinateNames(),'OutputName',outputFrameNames);
[estimated_parameters, estimated_delay] = parameterEstimation(r,data,parameterEstimationOptions);

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
