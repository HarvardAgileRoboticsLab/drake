function [xtraj, jltraj]  = variationalFourBar(plant,N,x0)

if nargin<1
% options.terrain = RigidBodyFlatTerrain();
% options.twoD = true; 
options.floating = false;
options.use_bullet = false;
%plant = PlanarRigidBodyManipulator(fullfile(getDrakePath,'matlab','systems','plants','test','FallingBrickContactPoints.urdf'),options);
%x0 = [0; 1; .1; 0; 0; 0];
file = fullfile(getDrakePath,'examples', 'SimpleFourBar', 'FourBar_JointLimits.urdf');
plant = RigidBodyManipulator(file,options);

% v = constructVisualizer(plant); 
% v.inspector([pi/2*ones(3,1); zeros(3,1)]); 

% add joint limits
% ub = [3; pi/6; pi/6];
% lb = [0; -pi/6; -pi/6]; 
% plant = plant.setJointLimits(lb, ub); 
% plant = compile(plant); 

end
if nargin < 2
  N=11;
end
if nargin<3
% q0 = [pi/6
%       0
%       0];
% q0 = [0; pi/6];
q0 = pi/2*ones(3,1);
% 
x0 = [q0;0*q0];
end

tf=2;
dt = tf/(N-1);

% ts_plant = TimeSteppingRigidBodyManipulator(plant,0.001,options);
sim_traj = plant.simulate([0,tf], x0);
% x0 = sim_traj.eval(0);

t_init = linspace(0,tf,N);
% traj_init.x = PPTrajectory(foh(sim_traj.getBreaks(), sim_traj.eval(sim_traj.getBreaks())));
traj_init.x = PPTrajectory(foh([0, tf], [x0, x0])); 
% 
options.s_weight = 10;
nq = plant.getNumPositions;
% 
traj_opt = VariationalTrajectoryOptimization(plant,N,tf,options);
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(x0(1:nq)),1);  
traj_opt = traj_opt.addVelocityConstraint(ConstantConstraint(x0(nq+(1:nq))),1);
% 
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',1000000);
traj_opt = traj_opt.setSolverOptions('snopt','MajorOptimalityTolerance',1e-7);
traj_opt = traj_opt.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-7);
traj_opt = traj_opt.setSolverOptions('snopt','MinorFeasibilityTolerance',1e-7);
traj_opt = traj_opt.setSolverOptions('snopt','constraint_err_tol',1e-5);

% traj_opt = traj_opt.setSolver('ipopt');

tic
[xtraj,~,~,~,~,~,jltraj, kltraj] = traj_opt.solveTraj(t_init,traj_init);
toc

% traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',1000000);
% traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',10000);
% traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
% traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',1000000);
% traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',1000);

v = constructVisualizer(plant);
v.playback(xtraj, struct('slider', true));
% 
% 
ts = sim_traj.getBreaks();
tt = xtraj.getBreaks();
h = tf/(N-1); 
xx_knot = xtraj.eval(tt); qq = xx_knot(1:nq, :); 
xx_mid = xtraj.eval(tt + h/2); vv = xx_mid(nq+1:end, :); 
xx = [qq; vv]; 
xs = sim_traj.eval(ts);

figure(1); clf; 
for i=1:size(xx, 1)
  subplot(2,size(xx,1)/2,i);
  if i <= nq
    plot(tt,xx(i,:),'b');
  else
      plot(tt + h/2, xx(i,:), 'b'); 
  end
  hold on;
  plot(ts,xs(i,:),'r--');
  hold off;
end

end

