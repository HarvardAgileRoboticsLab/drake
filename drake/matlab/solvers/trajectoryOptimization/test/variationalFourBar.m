function [xtraj,utraj,ctraj,btraj,psitraj,etatraj,jltraj, kltraj, straj]  ...
    = variationalFourBar(plant,N,x0,tf,x_init)

if nargin<1
    options.floating = false;
    options.use_bullet = false;
    file = fullfile(getDrakePath,'examples', 'SimpleFourBar4DOF', 'FourBar_JointLimits.urdf');
    plant = RigidBodyManipulator(file,options);    
end
if nargin < 2
    N=21;
end
if nargin<3
    q0 = [pi/4; 3*pi/4; pi/4]; %pi/4*ones(3,1);
    valuecheck(positionConstraints(plant, q0), zeros(6,1), 1e-4);
    x0 = [q0;0*q0];
end
if nargin<4
    tf = 2;
end
if nargin<5
    traj_init.x = PPTrajectory(foh([0 tf],[x0, x0]));
else
    traj_init.x = x_init; 
end

t_init = linspace(0,tf,N);


% ts_plant = TimeSteppingRigidBodyManipulator(plant,dt/2,options);
sim_traj = plant.simulate([0,tf], x0);
% x0 = sim_traj.eval(0);

options.s_weight = 10;
nq = plant.getNumPositions;

traj_opt = VariationalTrajectoryOptimization(plant,N,tf,options);
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(x0(1:nq)),1);
traj_opt = traj_opt.addVelocityConstraint(ConstantConstraint(x0(nq+(1:nq))),1);

% traj_opt = traj_opt.setSolver('ipopt');
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',1000000);
traj_opt = traj_opt.setSolverOptions('snopt','MajorOptimalityTolerance',1e-5);
traj_opt = traj_opt.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-5);
traj_opt = traj_opt.setSolverOptions('snopt','MinorFeasibilityTolerance',1e-5);
traj_opt = traj_opt.setSolverOptions('snopt','constraint_err_tol',1e-5);


tic
[xtraj,utraj,ctraj,btraj,psitraj,etatraj,jltraj, kltraj, straj] = traj_opt.solveTraj(t_init,traj_init);
toc

% traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',1000000);
% traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',10000);
% traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
% traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',1000000);
% traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',1000);

% v = constructVisualizer(plant);
% v.playback(xtraj, struct('slider', true));

dt = tf/N; 
ts = sim_traj.getBreaks();
tt = xtraj.getBreaks();
xx_knot = xtraj.eval(tt); qq = xx_knot(1:nq, :);
xx_mid = xtraj.eval(tt + dt/2); vv = xx_mid(nq+1:end, :);
xx = [qq; vv];
xs = sim_traj.eval(ts);

rms_err = rms(xx(:,1:end-1) - [interp1(ts', xs(1:nq,:)', tt(1:end-1)'), ...
    interp1(ts', xs(nq+1:2*nq,:)', tt(1:end-1)+dt/2', 'linear', 'extrap')]',2);
disp(rms_err);
% save(['rms_err_', num2str(N)], 'dt', 'rms_err');

figure(2); clf;
for i=1:size(xs, 1)
  subplot(2,size(xs,1)/2,i);
  if i <= nq
    plot(tt,rad2deg(xx(i,:)),'b');
  else
      plot(tt + dt/2, rad2deg(xx(i,:)), 'b');
  end
  hold on;
  plot(ts,rad2deg(xs(i,:)),'r--');
  hold off;
  legend('TrajOpt', 'TimeStepping')
end

end

