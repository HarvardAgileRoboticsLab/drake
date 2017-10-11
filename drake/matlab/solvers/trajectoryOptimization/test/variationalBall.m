function xtraj = variationalBall(plant,N,x0)

if nargin<1
options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.use_bullet = true;
%plant = PlanarRigidBodyManipulator(fullfile(getDrakePath,'matlab','systems','plants','test','FallingBrickContactPoints.urdf'),options);
%x0 = [0; 1; .1; 0; 0; 0];
file = fullfile(getDrakePath,'matlab','systems','plants','test','ball.urdf');
plant = RigidBodyManipulator(file,options);
end
if nargin < 2
  N=11;
end
if nargin<3
q0 = [0
      0
      1.0000
      0.1
     -0.1
      0.0];
    
x0 = [q0;0*q0];
end

tf=2.0;
dt = tf/(N-1);

%p = TimeSteppingRigidBodyManipulator(s,.1,options);
p = VariationalRigidBodyManipulator(plant,dt);

tic
xtraj = p.simulate([0 tf],x0);
toc

% v = p.constructVisualizer();
% v.playback(xtraj,struct('slider',true));

% t_init = linspace(0,tf,N);
% traj_init.x = PPTrajectory(foh([0 tf],[x0, x0]));
% 
% options.s_weight = 10;
% nq = plant.getNumPositions;
% 
% traj_opt = VariationalTrajectoryOptimization(plant,N,tf,options);
% traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(x0(1:nq)),1);  
% traj_opt = traj_opt.addVelocityConstraint(ConstantConstraint(x0(nq+(1:nq))),1);
% 
% %traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',1000000);
% traj_opt = traj_opt.setSolver('ipopt');
% 
% tic
% xtraj = traj_opt.solveTraj(t_init,traj_init);
% toc

traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',1000000);
traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',10000);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
% traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',1000000);
% traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',1000);
% traj_opt = traj_opt.setSolverOptions('snopt','MajorOptimalityTolerance',1e-7);
% % traj_opt = traj_opt.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-7);
% traj_opt = traj_opt.setSolverOptions('snopt','MinorFeasibilityTolerance',1e-7);
% traj_opt = traj_opt.setSolverOptions('snopt','constraint_err_tol',1e-5);


v = constructVisualizer(plant);
v.playback(xtraj);

% ts_plant = TimeSteppingRigidBodyManipulator(plant,0.0005,options);

% sim_traj = ts_plant.simulate([0,tf],x0);
% 
% 
% ts = sim_traj.getBreaks();
% %tt = xtraj.getBreaks();
% xx = xtraj.eval(ts);
% xs = sim_traj.eval(ts);
% 
% figure(1);
% for i=1:12
%   subplot(2,6,i);
%   plot(ts,xx(i,:),'b');
%   hold on;
%   plot(ts,xs(i,:),'r');
%   hold off;
% end

end

