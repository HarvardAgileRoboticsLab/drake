function [p,xtraj,utraj,ctraj,btraj,psitraj,etatraj,straj,z, ...
    F,info,infeasible_constraint_name, traj_opt] = runStepUp()

l = 0.4;
h = 0.11;
boxes = [0.25+l, 0.0, 2*l, 1, h];
%          0.25+l+l/2, 0.0, l, 1, 2*h];
options.terrain = RigidBodyStepTerrain(boxes);

options.floating = true;
options.ignore_self_collisions = true;
options.use_bullet = false;
p = LittleDog(options);
v = p.constructVisualizer;

nq = p.getNumPositions();
nv = p.getNumVelocities();
nx = nq+nv;
nu = p.getNumInputs();

% Load nominal data
x0 = double(home(p));
q0 = x0(1:nq);
q0(3) = q0(3) - 0.010;

T0 = 4;
N = 30;

% ----- Initial Guess ----- %
q1 = [0.45;0;q0(3)+h;q0(4:end)];
x1 = [q1;zeros(nv,1)];
q1_ub = q1 + 0.01;
q1_lb = q1 - 0.01;

% good_random_init = load('TimingResults/SFTimedTrial_good_ic.mat'); 
% 
t_init = linspace(0, T0, N); 
% t_init = [0; cumsum(good_random_init.t0)];
% traj_init.x = PPTrajectory(foh(t_init, [good_random_init.x0; 0*good_random_init.x0]));
% traj_init.u = PPTrajectory(zoh(t_init(1:end-1), good_random_init.u0));
% traj_init.c = PPTrajectory(zoh(t_init(1:end-1), good_random_init.c0));
% traj_init.b = PPTrajectory(zoh(t_init(1:end-1), good_random_init.b0));
% traj_init.eta = PPTrajectory(zoh(t_init(1:end-1), good_random_init.eta0));
% traj_init.psi = PPTrajectory(zoh(t_init(1:end-1), good_random_init.psi0));
% traj_init.s = PPTrajectory(zoh(t_init(1:end-1), good_random_init.s0'));

% traj_init.x = PPTrajectory(foh([0 T0/2 T0],[x0, xm, x1]));
traj_init.x = PPTrajectory(foh([0 T0],[x0, x1]));
traj_init.u = PPTrajectory(zoh(t_init,0.1*randn(nu,N)));
traj_init.s = PPTrajectory(zoh([0 T0], [1, 1])); 
T_span = [1 T0];

traj_opt = VariationalTrajectoryOptimization(p,N,T_span, struct('s_weight', 30));
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q0),1);  
% traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(qm),7);
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(q1_lb, q1_ub),N);
traj_opt = traj_opt.addVelocityConstraint(ConstantConstraint(zeros(nv,1)),1);

[q_lb, q_ub] = getJointLimits(p);
z_ub = q1(3) + 0.01; 
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(0, z_ub), 2:N-1, 3); 
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(q_lb,q_ub),2:N-1);


state_cost = Point(getStateFrame(p),ones(nx,1));
state_cost.base_x = 0;
state_cost.base_y = 0;
state_cost.base_z = 0;
state_cost.base_pitch = 10;
state_cost.base_roll = 10;
state_cost.base_yaw = 10;
state_cost.front_left_hip_roll = 5;
state_cost.front_right_hip_roll = 5;
state_cost.back_left_hip_roll = 5;
state_cost.back_right_hip_roll = 5;
% state_cost.front_left_hip_pitch = 5;
% state_cost.front_right_hip_pitch = 5;
% state_cost.back_left_hip_pitch = 5;
% state_cost.back_right_hip_pitch = 5;
state_cost = double(state_cost);
Q = diag(state_cost); 


% traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(qf(1:6)),N);
%traj_opt = traj_opt.addInputConstraint(ConstantConstraint(zeros(3,1)),1:N-1);
%traj_opt = traj_opt.addStateConstraint(ConstantConstraint(qm(1:6)),8);
% traj_opt = traj_opt.addStateConstraint(ConstantConstraint(qf(1:6)),N);
% traj_opt = traj_opt.addPositionConstraint(periodic_constraint,{[1 N]});
traj_opt = traj_opt.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-5);
traj_opt = traj_opt.setSolverOptions('snopt','MinorFeasibilityTolerance',1e-5);
traj_opt = traj_opt.setSolverOptions('snopt','MajorOptimalityTolerance',1e-5);
traj_opt = traj_opt.setSolverOptions('snopt','MinorOptimalityTolerance',1e-5);

traj_opt = traj_opt.setSolver('snopt');
traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',10000);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',1000000);
traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',1000);

% traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);

tic
[xtraj,utraj,ctraj,btraj,psitraj,etatraj,straj,z ...
    ,F,info,infeasible_constraint_name] = traj_opt.solveTraj(t_init,traj_init);
toc

v.playback(xtraj,struct('slider',true));

function [f,df] = running_cost_fun(h,x,u)
  R = 10*eye(nu);
  g = (1/2)*(x-x1)'*Q*(x-x1) + (1/2)*u'*R*u;
  f = h*g;
  df = [g, h*(x-x1)'*Q, h*u'*R];
end



  function displayTraj(h,x,u)
  
    ts = [0;cumsum(h)];
    for i=1:length(ts)
      v.drawWrapper(0,x(:,i));
      pause(h(1)/3);
    end
   
end
end
