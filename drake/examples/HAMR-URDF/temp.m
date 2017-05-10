clear; clc; close all; 

%% Build Actuator 
dp.Vb = 225; 
dp.Vg = 0; 

nact = 2; 
sl_actuators = HamrActuators(nact, {'FLlact', 'FLsact'}, [1; 1], dp);

%% Simulate Actuator Bank

dt = 0.001;     % sampling time
fd = 1;         % drive frequency (Hz)
tsim = 5/fd;    % simulation time

t = 0:dt:tsim;

Vact = [0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t + 3*pi/2);
        0.5*(dp.Vb-dp.Vg)*sin(2*pi*fd*t)]; 

% ramp
tramp = 3/fd;
ramp = t/tramp; ramp(t >= tramp) = 1;
Vact = bsxfun(@times, ramp, Vact) + 0.5*(dp.Vb - dp.Vg); 

u = PPTrajectory(foh(t, Vact));

u = u.setOutputFrame(sl_actuators.getInputFrame().getFrameByName('DriveVoltage'));

connection(1).from_output = u.getOutputFrame();
% connection(1).to_input = sl_actuators.getInputFrame(); 
connection(1).to_input = sl_actuators.getInputFrame().getFrameByName('DriveVoltage');

output_select(1).system = 2;
output_select(1).output = sl_actuators.getOutputFrame(); 


% connection1(2).from_output = u.getOutputFrame().getFrameByName('ActuatorDeflection');
% connection1(1).to_input = sl_actuators.getInputFrame(); 
% connection1(2).to_input = sl_actuators.getInputFrame().getFrameByName('ActuatorDeflection');

% r = cascade(u, u2); 
r = mimoCascade(u, sl_actuators, connection, [], output_select); 
ytraj = r.simulate([0, tsim]); 
Fact = ytraj.eval(t);

figure(1); clf; hold on;
% plot(t, Fact)
% yyaxis 
yyaxis left;  plot(t, Vact(1,:), 'b', t, Vact(2,:), 'b--'); 
yyaxis right;  plot(t, Fact(1,:), 'r', t, Fact(2,:), 'r--'); 
