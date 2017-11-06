clear; clc; close all;
load('sysid_traj.mat')
%% Build Full Single Leg

% options
name = 'FL_scaled';
urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'dev', 'SimpleHAMR', 'urdf', [name, '.urdf']);
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = false;
options.z_inactive_guess_tol = 0.1;
options.dt = mean(diff(t))/10;

SL = SLTSRBM(urdf, options);
nq = SL.getNumPositions();
nv = SL.getNumVelocities();
nu = SL.getNumInputs();

% Relax joint limits
SL = SL.setJointLimits(-Inf(nq, 1), Inf(nq, 1));
SL = compile(SL);

%% Build Simple Single Leg

sl_urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'dev', 'SimpleHAMR', ...
    'urdf','SLSimple_scaled.urdf');

% Initial guess for stiffness and damping
sfb_inputs = [3; 7];
xout = xi(sfb_inputs, :);
xin = xi(SL.getActuatedJoints(), :);
Trat = xin'\xout';
Trat_sq = Trat*Trat;

K0 = zeros(numel(sfb_inputs)); %-diag(3*0.05*ones(2,1)) - Trat_sq\diag([1.1556; 1.1556]);
P0 = zeros(numel(sfb_inputs)); %-diag(0.01*ones(2,1)) - Trat_sq\diag([0.01; 0.01]);

options.k = K0(:);
options.p = P0(:);

SLSimple = SLSimpleRBM(sl_urdf, options);
viz = SLSimple.constructVisualizer();

x0 = SLSimple.getInitialState();
nqS = SLSimple.getNumPositions();
nvS = SLSimple.getNumVelocities();

tau_traj = PPTrajectory(foh(t, tau));
tau_traj = tau_traj.setOutputFrame(SLSimple.getInputFrame());

% Simulate with initial guess
i_start = find(t >= 1000, 1, 'first');
i_end = i_start+1000;

SLSimple_OL0 = cascade(tau_traj, SLSimple);
xtraj0 = simulate(SLSimple_OL0, [t(1), t(i_end)], xi([sfb_inputs; nq+sfb_inputs], 1));
xx0 = xtraj0.eval(t);

%% Fitting

params.SL = SL;
params.SLSimple = SLSimple;
params.t = t(i_start:i_end);
params.x = xi(1:(nq+nv), i_start:i_end);
params.u = u(:, i_start:i_end);
params.tau = tau(:, i_start:i_end);
params.sfb_inputs = sfb_inputs;

z0 = [K0(:); P0(:)];
lb = -Inf(size(z0)); %[0; -Inf; -Inf; 0; 0; -Inf; -Inf; 0];
ub = Inf(size(z0));

optopt = optimoptions('lsqnonlin', 'SpecifyObjectiveGradient', true);
[zopt,resnorm,residual,exitflag,output] = lsqnonlin(@(z0)SLFitFun(z0, params),z0,lb,ub,optopt);

SLSimple = SLSimple.setK(reshape(zopt(1:nqS^2), nqS, nqS)); disp(SLSimple.K);
SLSimple = SLSimple.setP(reshape(zopt(nqS^2+(1:nvS^2)), nvS, nvS)); disp(SLSimple.P);
SLSimple = SLSimple.removeAllStateConstraints();

%% Verify
SLSimple_OL = cascade(tau_traj, SLSimple);
xtraj_fit = simulate(SLSimple_OL, [t(1), t(i_end)], xi([sfb_inputs; nq+sfb_inputs], 1));
xxfit = xtraj_fit.eval(t);
%% Plotting
ind_plot = [sfb_inputs; nq+sfb_inputs];
figure(1); clf; hold on;
titles= {'Swing', 'Lift', 'Swing Rate', 'Lift Rate'};
for i = 1:numel(ind_plot)
    subplot(2,2,i); hold on;
    title(titles{i});
    plot(t, xi(ind_plot(i), :))
    plot(t, xx0(i, :) - x0(i));
    plot(t, xxfit(i, :) - x0(i))
end

xtraj_fit_scaled = PPTrajectory(foh(t*1e-3, xxfit));
xtraj_fit_scaled = xtraj_fit_scaled.setOutputFrame(xtraj_fit.getOutputFrame());
viz.playback(xtraj_fit_scaled, struct('slider', true));

rms(xi(ind_plot, :) - (xxfit - x0),2)