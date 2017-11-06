clear; clc; close all;

% Sample Rate and Time
dt = 1;
Fs = 1/dt;
T = 5000;
t = 0:dt:T;
N = numel(t);

% White Noise Freq Characteristics
f1 = 150e-3;
band = [0, 2*f1*dt];

%% Build robot

% options
name = 'FL_scaled';
urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'dev', 'SimpleHAMR', 'urdf', [name, '.urdf']);
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = false;
options.z_inactive_guess_tol = 0.1;
options.dt = dt/10;

SL = SLTSRBM(urdf, options);
nq = SL.getNumPositions();
nv = SL.getNumVelocities();
nu = SL.getNumInputs();

% Relax joint limits
SL = SL.setJointLimits(-Inf(nq, 1), Inf(nq, 1));
SL = compile(SL); 


%% Build Actuators
dp.Vb = 100;
dp.Vg = 0;

nact = nu;
sl_actuators = HamrActuators(nact, {'FLsact', 'FLlact'}, [1; 1], dp);

%% Connect system

% connections from actuators to SL
sl_actuators = sl_actuators.setOutputFrame(SL.getInputFrame());
connection1(1).from_output = sl_actuators.getOutputFrame();
connection1(1).to_input = SL.getInputFrame();

% connections from SL to actuators
sl_out = SL.getOutputFrame();
act_in = sl_actuators.getInputFrame();
act_in = act_in.replaceFrameNum(2, sl_out.getFrameByName('ActuatorDeflection'));
sl_actuators = sl_actuators.setInputFrame(act_in);
%
connection2(1).from_output = sl_out.getFrameByName('ActuatorDeflection');
connection2(1).to_input = act_in.getFrameByName('ActuatorDeflection');

% mimo inputs
input_select(1).system = 1;
input_select(1).input = act_in.getFrameByName('DriveVoltage');

% mimo outputs
output_select(1).system = 2;
output_select(1).output = sl_out.getFrameByName('FrontLeftSingleLegPosition');
output_select(2).system = 2;
output_select(2).output = sl_out.getFrameByName('FrontLeftSingleLegVelocity');
output_select(3).system = 1;
output_select(3).output = sl_actuators.getOutputFrame();

SLWact = mimoFeedback(sl_actuators, SL, connection1, connection2, ...
    input_select, output_select);

%% Build and Plot Input

tramp = 1000; 
Vscale = (dp.Vb - dp.Vg)/5;
V = Vscale*idinput([N, SL.getNumInputs()], 'rgs', band)'; 
V(t < 1000) = (t(t<1000)/tramp).*V(t<1000);
V = V + (dp.Vb - dp.Vg)/2 + dp.Vg;
V(V > dp.Vb) = dp.Vb;
V(V < dp.Vg) = dp.Vg;

figure(1); clf;
subplot(2,1,1);
yyaxis left; hold on; plot(t, V(1,:))
subplot(2,1,2);
yyaxis left;  hold on; plot(t, V(2,:))

freq = [-(Fs/2):(1/N/dt):(Fs/2)];

f0 = 0;
ind0 = find(freq > f0, 1, 'first');
ind1 = find(freq > 2*f1, 1, 'first');
Vfft = abs(fftshift((fft(bsxfun(@minus, V, (dp.Vb - dp.Vg)/2 + dp.Vg), [], 2)/size(V,2))));

figure(2); clf;
subplot(2,1,1); 
yyaxis left; hold on; plot(freq(ind0:ind1)*1e3, 2*Vfft(1,ind0:ind1))
subplot(2,1,2); hold on;
yyaxis left; hold on; plot(freq(ind0:ind1)*1e3, 2*Vfft(2,ind0:ind1))

%% IC and  Simulate

Vtraj = PPTrajectory(zoh(t, V));
Vtraj = Vtraj.setOutputFrame(SLWact.getInputFrame);
SLWact_OL = cascade(Vtraj, SLWact);

q0 = zeros(nq,1);
x0 = [q0; 0*q0];

disp('Simulating...'); 
tic; xtraj = simulate(SLWact_OL, [0 T], x0); tlcp = toc;

tsim = xtraj.getBreaks();
x = xtraj.eval(tsim);

%%  Plot Response
sfb_inputs = [3; 7];
y = x(sfb_inputs, :);
yi = interp1(tsim', y', t')';
yfft = abs(fftshift((fft(yi, [], 2)/size(y,2))));

figure(1);
subplot(2,1,1); 
yyaxis right; hold on; plot(tsim, rad2deg(y(1,:)));
yyaxis right; hold on; plot(t, rad2deg(yi(1,:)));
subplot(2,1,2); hold on;
yyaxis right; hold on; plot(tsim, rad2deg(y(2,:)));
yyaxis right; hold on; plot(t, rad2deg(yi(2,:)));

figure(2);
subplot(2,1,1);  hold on;
yyaxis right; hold on;plot(freq(ind0:ind1)*1e3, 2*yfft(1,ind0:ind1))
subplot(2,1,2); hold on;
yyaxis right; hold on;plot(freq(ind0:ind1)*1e3, 2*yfft(2,ind0:ind1))

%% Convert Voltage to Torque about Leg

xi = interp1(tsim', x', t')';
u = zeros(nu, N); 
tau = zeros(nu, N);

actuated_dof = SL.getActuatedJoints();
for i = 1:N
    q = xi(1:nq, i);
    qd = xi(nq+(1:nv), i); 
    kinsol = SL.doKinematics(q, qd, struct('compute_gradients', true));
    [K, dK] = SL.positionConstraints(q);
    dKUnique = dK([1:2; 8:9; 13, 15], :);
    Jc = -dKUnique(:,[2:4, 6:8])\dKUnique(:, actuated_dof);
    Jc_phi = Jc([2, 5], :);
    tau(:,i) = Jc_phi\xi(nq+nv+(1:nu),i);    
end

figure(3); clf;
subplot(2,1,1); hold on; 
plot(t, xi(nq+nv+1,:)); 
plot(t, tau(1,:)); 

subplot(2,1,2); hold on; 
plot(t, xi(nq+nv+2,:)); 
plot(t, tau(2,:)); 

save('sysid_traj', 't', 'xi', 'u', 'tau')
