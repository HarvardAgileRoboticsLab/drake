%SIMULATIONINITFCN Simulation model initialization callback function.
%
% Syntax:
%   simulationInitFcn;

% Copyright 2015-2017 Mikhail S. Jones

% Clean up workspace 
clear all; close all; clc;

% Start clock at initialization (used for computing real time performance)
initTime = tic;

% Initialize data structures and interface objects
simulationParameters = CassieSimulationParameters;
cassieParameters = CassieParameters;
cassieEtherCat = CassieEtherCat;
cassieInputs = CassieInputs;
cassieOutputs = CassieOutputs;
cassieCalibration = CassieCalibration;

% Construct inverse kinematics problem
problem = NonlinearLeastSquaresProblem(...
  InverseKinematicsFunction, ...
  CassieParameters.motorPositionLowerLimit, ...
  CassieParameters.motorPositionUpperLimit);

% Construct inverse kinematics solver
solver = NonlinearLeastSquaresSolver(problem);

% Set desired foot frames
problem.residualFunction.setDesiredFootTransforms(...
  Transform3d().translate([0; 0.135; -1]), ...
  Transform3d().translate([0; -0.135; -1]));

% Solve inverse kinematics
solver.solve;

% Set initial motor positions
cassieEtherCat.setMotorPositions(solver.getSolution);

% Set initial pelvis position and rotation
pelvisPosition = [0; 0; 1.01];
pelvisRotation = Rotation3d().rotZYX([0; 0; 0]).getValue;

% Initialize Cassie EtherCAT bus signals
cassieEtherCatData = cassieEtherCat.getStructure;
cassieEtherCatBusInfo = Simulink.Bus.createObject(cassieEtherCatData);
cassieEtherCatBus = eval(cassieEtherCatBusInfo.busName);

% Initialize Cassie inputs bus signals
cassieInputsData = cassieInputs.getStructure;
cassieInputsBusInfo = Simulink.Bus.createObject(cassieInputsData);
cassieInputsBus = eval(cassieInputsBusInfo.busName);

% Initialize Cassie outputs bus signals
cassieOutputsData = cassieOutputs.getStructure;
cassieOutputsBusInfo = Simulink.Bus.createObject(cassieOutputsData);
cassieOutputsBus = eval(cassieOutputsBusInfo.busName);

% Initialize blank cassie calibration file that lives on target PC
cassieCalibration.saveData;