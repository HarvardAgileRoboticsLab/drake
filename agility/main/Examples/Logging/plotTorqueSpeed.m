%PLOTTORQUESPEED Plot torque speed curves from logged data.

% Copyright 2017 Agility Robotics
% Author: Mikhail Jones

% Run script to load log files
loadFileScopeData;

% Define motor specs
maxMotorVelocity = CassieParameters.motorMaximumSpeed./CassieParameters.motorGearRatio*2*pi/60;
maxMotorTorque = CassieParameters.motorMaximumTorque.*CassieParameters.motorGearRatio;

% Resample motor torques and positions
t = outputsTime;
dq = motorVelocity;
u = interp1(inputsTime, motorTorque, t);

% Plot
figure;
for i = 1:5
  subplot(2,3,i); hold on; grid on; box on;
  xlabel('Motor Output Velocity (rad/s)');
  ylabel('Motor Output Torque (N*m)');
  plot(dq(:,i), u(:,i));
  plot(dq(:,i + 5), u(:,i + 5));
  plot([-maxMotorVelocity(i) 0 maxMotorVelocity(i) 0 -maxMotorVelocity(i)], [0 maxMotorTorque(i) 0 -maxMotorTorque(i) 0], '--r');
end % for