clear; clc; close all;


dataset = '_80Hz_100V';
data = load(['sysid_traj', dataset, '.mat']);

N0 = find(data.t>=500, 1, 'first');
N = 1000;
ds = (numel(data.t) - N0)/N; 
IND = N0:ds:numel(data.t) ;      % randomly choose indices
SFB_INPUTS = [3;7];

NK = 2; 
NP = 1; 

[xf,objval,exitflag,infeasible_constraint_name] = SLSimpleFitNLP(data, IND, SFB_INPUTS, NK, NP); 

%% Simulate w/ fitted params

sl_urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'dev', 'SimpleHAMR', ...
    'urdf','SLSimple_scaled.urdf');

options.k = reshape(xf(1:4*NK), [], NK);
options.p = reshape(xf(4*NK+(1:4*NP)), [], NP); 


SLSimple = SLSimpleRBM(sl_urdf, options);
v = SLSimple.constructVisualizer();
nqS = SLSimple.getNumPositions();
nvS = SLSimple.getNumVelocities();
x0 = SLSimple.getInitialState();

K = reshape(options.k, nqS, nqS, NK);
P = reshape(options.p, nvS, nvS, NP);
save(['SpringDamper', dataset '_', num2str(NK), '_', num2str(NP)], 'K', 'P', 'objval'); 

t = data.t; 
tau_traj = PPTrajectory(foh(t, data.tau));
tau_traj = tau_traj.setOutputFrame(SLSimple.getInputFrame());
SLSimple_OL0 = cascade(tau_traj, SLSimple);

xtraj = simulate(SLSimple_OL0, [t(1), t(end)], x0);

%% Validate 

xx = xtraj.eval(t);
xtraj_scaled = PPTrajectory(foh(t*1e-3, xx));
xtraj_scaled = xtraj_scaled.setOutputFrame(xtraj.getOutputFrame());
v.playback(xtraj_scaled, struct('slider', true));

xd = data.x;
ind_plot = [SFB_INPUTS; 8+SFB_INPUTS];
figure(1); clf; hold on;
titles= {'Swing', 'Lift', 'Swing Rate', 'Lift Rate'};
for i = 1:numel(ind_plot)
    subplot(2,2,i); hold on;
    title(titles{i});
    plot(t, rad2deg(xd(ind_plot(i), :)))
    plot(t, rad2deg(xx(i, :) - x0(i)));
end

save(['SpringDamper', dataset '_', num2str(NK), '_', num2str(NP)], 'K', 'P', 'xf', 'objval'); 


