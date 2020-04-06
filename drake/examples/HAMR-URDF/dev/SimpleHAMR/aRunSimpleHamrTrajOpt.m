clear; clc; close all;

save_dir = './traj_opts/';
trial_name = 'straight_non_periodic12';

tht_vec = deg2rad(0:30:360);

% parameters
save_path = './traj_opts/straight/';
load_path = './traj_opts/straight_non_periodic11';
params.r_des = 50;
params.tht_des = 0;
params.N = 41;
params.s_weight = 70;

[hamr,x,u,c,b,psi,eta,jl, kl, s, z,F,info,infeasible_constraint_name] = ...
    SimpleHAMRVariationalTrajOpt(load_path, params);

save([save_dir, trial_name], 'x', 'u', 'c', 'b', 'psi', 'eta', ...
    'jl', 'kl', 's', 'params')
%%
v = hamr.constructVisualizer();
tt = x.getBreaks();
xx = x.eval(tt);
cc = c.eval(tt);

nq = hamr.getNumPositions;
x_playback = PPTrajectory(foh(tt/1e3, xx(1:nq,:)));
x_playback = x_playback.setOutputFrame(v.getInputFrame);
v.playback(x_playback, struct('slider', true));

%%
figure(1); clf;
for i = 1:6
    subplot(2,3,i); hold on;
    if i <= 3
        plotyy(tt*1e-3, xx(i,:), tt*1e-3, xx(i+nq, :))
    else
        plotyy(tt*1e-3, rad2deg(xx(i,:)), tt*1e-3, rad2deg(xx(i+nq,:)))
    end
end

figure(2); clf;
subplot(2,1,1); hold on;
plot(tt*1e-3, cc(1,:), tt*1e-3, cc(4,:))
subplot(2,1,2); hold on;
plot(tt*1e-3, cc(2,:), tt*1e-3, cc(3,:))