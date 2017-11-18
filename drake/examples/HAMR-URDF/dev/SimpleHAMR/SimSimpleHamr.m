clear; clc; close all;

%% Simple HAMR
options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = true;
options.collision = true;
options.dt = 1;

urdf_simple = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'dev', 'SimpleHAMR', ...
    'urdf', 'HAMRSimple_scaled.urdf');

hamr_simpleTSRBM = HamrSimpleTSRBM(urdf_simple, options);
v = hamr_simpleTSRBM.constructVisualizer();

%% Drive Signal

V = 100;
V2Tau = 0.0013; % ratio of voltage to torque

fd = 0.010;         % drive frequency (Hz)
tsim = 0.5e3;

t = 0:options.dt:tsim;

tau_act =  V*V2Tau*[sin(2*pi*fd*t + pi/2);            % FLswing
    sin(2*pi*fd*t);                                   % FLlift
    sin(2*pi*fd*t - pi/2);                                   % RLSwing
    sin(2*pi*fd*t + pi);                              % RLLift
    sin(2*pi*fd*t - pi/2);                            % FRswing
    sin(2*pi*fd*t);                              % FRlift
    sin(2*pi*fd*t + pi/2);                            % RRSwing
    sin(2*pi*fd*t + pi)];                                  % RRLift

% ramp
tramp = 1/fd;
ramp = t/tramp; ramp(t >= tramp) = 1;

tau_act = bsxfun(@times, ramp, tau_act);
u = PPTrajectory(foh(t, tau_act));
u = setOutputFrame(u, hamr_simpleTSRBM.getInputFrame());

%% Simulate Open loop

hamr_OL = cascade(u, hamr_simpleTSRBM);
disp('Valid initial condition: simulating...')
xtraj = simulate(hamr_OL, [0 tsim], hamr_simpleTSRBM.getInitialState());
xtraj_scaled = DTTrajectory(xtraj.getBreaks()*1e-3, xtraj.eval(xtraj.getBreaks()));
xtraj_scaled = xtraj_scaled.setOutputFrame(xtraj.getOutputFrame());

%%

tt = xtraj_scaled.getBreaks();
xx = xtraj_scaled.eval(tt);
nq = hamr_simpleTSRBM.getNumPositions();

state_frame = Point(hamr_simpleTSRBM.getStateFrame());

if options.floating
    i0 = 7;
else
    i0 = 1;
end
for i = i0:nq
    figure(i-i0+1); clf; hold on;    
    plot(tt, rad2deg(xx(i,:)))
end
tilefigs;

v = hamr_simpleTSRBM.constructVisualizer();
options.slider = true;
v.playback(xtraj_scaled, options);


