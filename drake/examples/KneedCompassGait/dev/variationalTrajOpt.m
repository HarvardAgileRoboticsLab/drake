function [p,xtraj,utraj,z,F,info,traj_opt] = variationalTrajOpt()

%[~,xtraj,utraj] = testNewTrajOpt();

options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.ignore_self_collisions = true;
options.use_bullet = false;
p = PlanarRigidBodyManipulator('../KneedCompassGait.urdf',options);

N = 17;
T0 = 3;

% ----- Periodicity Constraint ----- %
% Positions
q_periodic = zeros(p.getNumPositions);
q_periodic(2,2) = 1; %z
q_periodic(3,3) = 1; %pitch-hip w/symmetry
q_periodic(3,5) = 1; %pitch-hip w/symmetry
q_periodic(4,6) = 1; %knee w/symmetry
q_periodic(6,4) = 1; %knee w/symmetry
q_periodic(5,5) = -1; %hip w/symmetry

% %Velocities
% v_periodic = zeros(p.getNumPositions);
% v_periodic(1,1) = 1; %x-vel
% v_periodic(2,2) = 1; %z-vel
% v_periodic(3,3) = 1; %pitch-hip w/symmetry
% v_periodic(3,5) = 1; %pitch-hip w/symmetry
% v_periodic(4,6) = 1; %knee w/symmetry
% v_periodic(6,4) = 1; %knee w/symmetry
% v_periodic(5,5) = -1; %hip w/symmetry
% 
% %Difference matrix to calculate velocities
% Diff = [-eye(p.getNumPositions), eye(p.getNumPositions)];
% 
% %Constraint Matrix
% R_periodic = [q_periodic, zeros(p.getNumPositions), zeros(p.getNumPositions), blkdiag(0, -eye(p.getNumPositions-1));
%               v_periodic*Diff, -eye(p.getNumPositions)*Diff];
          
R_periodic = [q_periodic, blkdiag(0, -eye(p.getNumPositions-1))];
              
periodic_constraint = LinearConstraint(zeros(p.getNumPositions,1),zeros(p.getNumPositions,1),R_periodic);
periodic_constraint = periodic_constraint.setName('periodicity');

% ----- Initial Guess ----- %

% x0 = [0;1;zeros(10,1)];
% x1 = [.1;1;pi/8-pi/16;pi/8;-pi/8;pi/8;zeros(6,1)];
% xf = [.2;1;zeros(10,1)];

q0 = [-.2 1 0 0 0 0]';
qm = [0 1 -.2 .2 .4 -.2]';
q1 = [.2 1 0 0 0 0]';
x0 = [q0; zeros(6,1)];
xm = [qm; zeros(6,1)];
x1 = [q1; zeros(6,1)];


t_init = linspace(0,T0,N);
% traj_init.x = PPTrajectory(foh([0 T0/2 T0],[x0 x1 xf]));
traj_init.x = PPTrajectory(foh([0 T0/2 T0],[x0, xm, x1]));
traj_init.u = PPTrajectory(zoh(t_init,0.1*randn(3,N)));
T_span = [1 T0];

traj_opt = VariationalTrajectoryOptimization(p,N,T_span);
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q0(1:6)),1);  
% traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(qm(1:6)),7);
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q1(1:6)),N);
traj_opt = traj_opt.addVelocityConstraint(ConstantConstraint(zeros(6,1)),1);

q_lb = [-inf; -inf; -inf; 0; -1.0; 0];
q_ub = [ inf;  inf;  inf; pi/2; 1.0; pi/2];
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(q_lb,q_ub),2:N-1);

% traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(qf(1:6)),N);
%traj_opt = traj_opt.addInputConstraint(ConstantConstraint(zeros(3,1)),1:N-1);
%traj_opt = traj_opt.addStateConstraint(ConstantConstraint(qm(1:6)),8);
% traj_opt = traj_opt.addStateConstraint(ConstantConstraint(qf(1:6)),N);
% traj_opt = traj_opt.addPositionConstraint(periodic_constraint,{[1 N]});

traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',10000);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',1000000);

tic
[xtraj,utraj,ctraj,btraj,straj,z,F,info,infeasible_constraint_name] = traj_opt.solveTraj(t_init,traj_init);
toc

v = p.constructVisualizer();
v.playback(xtraj,struct('slider',true));

function [f,df] = running_cost_fun(h,x,u)
  Q = 100*eye(12);
  R = 1*eye(3);
  g = (1/2)*(x-x1)'*Q*(x-x1) + (1/2)*u'*R*u;
  f = h*g;
  df = [g, h*(x-x1)'*Q, h*u'*R];
end

end
