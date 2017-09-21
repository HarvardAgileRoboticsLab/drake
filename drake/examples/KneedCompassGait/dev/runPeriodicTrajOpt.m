clear; clc; close all;

[plant, xtraj,utraj,ctraj,btraj, psitraj,etatraj,jltraj, kltraj,straj, ...
    z,F,info,infeasible_constraint_name] = variationalTrajOptPeriodic();

qtraj = xtraj(1:6);
qtraj = qtraj.setOutputFrame(v.getInputFrame());
v.playback(qtraj, struct('slider', 'true'))

tt = xtraj.getBreaks();
N = numel(tt);
T_span = [tt(1), tt(end)];
traj_opt = VariationalTrajectoryOptimization(plant,N,T_span);
h = mean(diff(tt)); 

xx = xtraj.eval(tt); 
uu = utraj.eval(tt); 
cc = ctraj.eval(tt); 
bb = btraj.eval(tt);
% jj = jltraj.eval(tt);
% kk = kltraj.eval(tt);

p0 = left_legendre_transform_fun(traj_opt,h,xx(1:6,1),xx(1:6,2),uu(:,1), ...
    cc(:,1),bb(:,1),[],[]);
pN = right_legendre_transform_fun(traj_opt,h,xx(1:6,end-1),xx(1:6,end),uu(:,end-1),[]);
dp = p0 - pN

M0 = plant.manipulatorDynamics(xx(1:6,1), zeros(6,1)); 
MN = plant.manipulatorDynamics(xx(1:6,end), zeros(6,1)); 

v0 = M0\p0; 
vN = MN\pN;
dv = vN - v0

