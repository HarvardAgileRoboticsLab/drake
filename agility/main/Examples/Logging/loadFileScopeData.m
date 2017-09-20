%LOADFILESCOPEDATA Load real-time file scope data into workspace.

% Copyright 2017 Agility Robotics
% Author: Mikhail Jones

% Message Scope -----------------------------------------------------------

% Display info message
fprintf('Loading file scope data from message.dat\n');

% Load data from realtime file scope
fileScopeData = SimulinkRealTime.utils.getFileScopeData('message.dat');

% Parse data
messageTime = fileScopeData.data(:,end);
message = fileScopeData.data(:,1:end-1);

% Initialize figure
figure;

% Plot message codes
messagePlotHandle = subplot(1,1,1);
hold on; grid on; xlabel('Time (s)'); ylabel('Message Code');
stem(messageTime, message);

% % Find all error messages
% [r, c] = find(message);
% 
% % Loop through found message
% for j = 1:length(r)
%   % Get current message properties
%   t = messageTime(r(j));
%   msg =  char(CassieState(message(r(j),c(j))));
% 
%   % Display time and error message
%   fprintf('(%0.4f) %s\n', t, msg);
% end % for


% Energy Scope ------------------------------------------------------------

% Display info message
fprintf('Loading file scope data from energy.dat\n');

% Load data from realtime file scope
fileScopeData = SimulinkRealTime.utils.getFileScopeData('energy.dat');

% Find the first nonzero value
start = find(any(fileScopeData.data(:,1:end-1), 2), 1);

% Parse data
energyTime = fileScopeData.data(start:end,end);
V = fileScopeData.data(start:end,1);
A = fileScopeData.data(start:end,2);

% Compute power draw
P = V.*A;

% Initialize figure
figure;

% Plot voltage
voltagePlotHandle = subplot(3, 1, 1);
hold on; grid on; xlabel('Time (s)'); ylabel('Voltage (V)');
plot(energyTime, V);

% Plot current
currentPlotHandle = subplot(3, 1, 2);
hold on; grid on; xlabel('Time (s)'); ylabel('Current (A)');
plot(energyTime, A);

% Plot power
powerPlotHandle = subplot(3, 1, 3);
hold on; grid on; xlabel('Time (s)'); ylabel('Power (W)');
plot(energyTime, P);

% Link subplot axes
linkaxes([voltagePlotHandle, currentPlotHandle, powerPlotHandle], 'x');


% Outputs Scope -----------------------------------------------------------

% Display info message
fprintf('Loading file scope data from outputs.dat\n');

% Load data from realtime file scope
fileScopeData = SimulinkRealTime.utils.getFileScopeData('outputs.dat');

% Find the first nonzero value
start = find(any(fileScopeData.data(:,1:end-1), 2), 1);

% Parse data
outputsTime = fileScopeData.data(start:end,end);
vectorNavOrientation = fileScopeData.data(start:end,1:4);
vectorNavVelocity = fileScopeData.data(start:end,5:7);
motorPosition = fileScopeData.data(start:end,8:17);
motorVelocity = fileScopeData.data(start:end,18:27);
jointPosition = fileScopeData.data(start:end,28:31);
jointVelocity = fileScopeData.data(start:end,32:35);

% Initialize figure
figure;

% Plot motor position
motorPositionPlotHandle = subplot(2,1,1);
hold on; grid on; xlabel('Time (s)'); ylabel('Motor Position (rad)');
plot(outputsTime, motorPosition);

% Plot motor velocity
motorVelocityPlotHandle = subplot(2,1,2);
hold on; grid on; xlabel('Time (s)'); ylabel('Motor Velocity (rad/s)');
plot(outputsTime, motorVelocity);

% Initialize figure
figure;

% Plot joint position
jointPositionPlotHandle = subplot(2,1,1);
hold on; grid on; xlabel('Time (s)'); ylabel('Joint Position (rad)');
plot(outputsTime, jointPosition);

% Plot joint velocity
jointVelocityPlotHandle = subplot(2,1,2);
hold on; grid on; xlabel('Time (s)'); ylabel('Joint Velocity (rad/s)');
plot(outputsTime, jointVelocity);


% Inputs Scope ------------------------------------------------------------

% Display info message
fprintf('Loading file scope data from inputs.dat\n');

% Load data from realtime file scope
fileScopeData = SimulinkRealTime.utils.getFileScopeData('inputs.dat');

% Find the first nonzero value
start = find(any(fileScopeData.data(:,1:end-1), 2), 1);

% Parse data
inputsTime = fileScopeData.data(start:end,end);
motorTorque = fileScopeData.data(start:end,1:10);

% Initialize figure
figure;

% Plot  data
motorTorquePlotHandle = subplot(1,1,1);
hold on; grid on; xlabel('Time (s)'); ylabel('Motor Torque (N*m)');
plot(inputsTime, motorTorque);


% User Interface ----------------------------------------------------------

% Link output and input plots
linkaxes([...
  motorPositionPlotHandle, ...
  motorVelocityPlotHandle, ...
  jointPositionPlotHandle, ...
  jointVelocityPlotHandle, ...
  motorTorquePlotHandle], 'x');