%SIMULATIONSTOPFCN Simulation model stop callback function.
%
% Syntax:
%   simulationStopFcn;

% Copyright 2015-2017 Mikhail S. Jones

% Store computation time
computationTime = toc(initTime) - buildTime;

% Get simulation time
simulationTime = get_param(bdroot, 'SimulationTime');

% Compute real time performance
percentRealTime = 100*simulationTime/computationTime;

% Display
fprintf('Build time: %0.2f seconds \n', buildTime);
fprintf('Simulation time: %0.2f seconds \n', simulationTime);
fprintf('Computation time: %0.2f seconds \n', computationTime);
fprintf('Performance: %0.2f%% real time \n', percentRealTime);