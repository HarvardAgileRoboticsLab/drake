clear; clc; close all; 

% options.twoD = true;
options.floating = true; 
options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;

dt = 5e-5; 
Fs = 1/dt; 
tsim = 2; 
t = 0:dt:tsim; 

p = TimeSteppingRigidBodyManipulator('vibrationfoot2d.urdf', dt, options);
v = p.constructVisualizer();
v.xlim = [-.5 .5]; v.ylim = [-0.1 0.9];

nq = p.getNumPositions();
% q0 = 0; 
q0 = [0; 0; 0.005; 0; 0; 0; 0];

%% Generate white noise input

fmax = 500;
band = [0, 2*500]*dt;

urms = 100; 
uu = urms*idinput(numel(t), 'rgs', band)'; 
u = PPTrajectory(foh(t, uu));
u = setOutputFrame(u, p.getInputFrame());

%% Simulate 

pff = cascade(u, p);
xtraj = pff.simulate([0 tsim], [q0; 0*q0]);
tt = xtraj.getBreaks();
xx = xtraj.eval(xtraj.getBreaks());

%% Temoporal

ns = 2;

figure(1); clf; 
s1 = subplot(3,1,1); hold on; 
set(s1, 'fontsize', 14);
plot(t, uu); ylabel('Motor Force (N)')

s2 = subplot(3,1,2); hold on; 
set(s2, 'fontsize', 14);
ylabel('Foot Height (cm)')
plot(tt(1:ns:end), 10*xx(3,1:ns:end));

s3 = subplot(3,1,3); hold on; 
set(s3, 'fontsize', 14);
plot(tt(1:ns:end), 10*xx(7,1:ns:end));
ylabel('Spring Length (cm)')
xlabel('Time (s)')

print('-dpdf', '-bestfit', ['TimeData_urms=', num2str(urms)])

%% Frequency
 
freq = -(Fs/2):(1/numel(t)/dt):(Fs/2);
xxfft = abs(fftshift((fft(xx, [], 2)/size(xx,2))));
ufft = abs(fftshift((fft(uu, [], 2)/size(uu,2))));

ind0 = find(freq >= 0, 1, 'first');
ind1 = find(freq >= fmax, 1, 'first');


figure(2); clf; 
s1 = subplot(3,1,1); hold on; 
set(s1, 'fontsize', 14);
plot(freq(ind0:ind1), 2*ufft(1,ind0:ind1)); ylabel('Motor Force (N)')
xlim([0, fmax])

s2 = subplot(3,1,2); hold on; 
set(s2, 'fontsize', 14);
plot(freq(ind0:ind1), 2*xxfft(3,ind0:ind1)); 
ylabel('Foot Height (cm)')
xlim([0, fmax])

s3 = subplot(3,1,3); hold on; 
set(s3, 'fontsize', 14);
plot(freq(ind0:ind1), 2*xxfft(7,ind0:ind1));
ylabel('Spring Length (cm)')
xlabel('Frequency (Hz)')
xlim([0, fmax])

print('-dpdf', '-bestfit', ['FreqData_urms=', num2str(urms)])

%% 
% tilefigs; 
v.playback(xtraj, struct('slider', true));

