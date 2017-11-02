function [p,v,xtraj,utraj,ctraj,btraj,psitraj,etatraj,straj] = runVariationalTrajOpt(x_guess,u_guess,c_guess,b_guess,psi_guess,eta_guess,s_guess)

options.terrain = RigidBodyFlatTerrain();
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

T0 = 2;
N = 10;

% ----- Initial Guess ----- %
q1 = [0.6;q0(2:end)];
x1 = [q1;zeros(nv,1)];

if nargin > 0
    traj_init.x = x_guess;
    traj_init.u = u_guess;
    traj_init.c = c_guess;
    traj_init.b = b_guess;
    if nargin > 4
        traj_init.psi = psi_guess;
        traj_init.eta = eta_guess;
    end
    if nargin > 6
        traj_init.s = s_guess;
    end
    t_init = x_guess.getBreaks();
    if length(t_init) ~= N
        t_init = linspace(0,t_init(end),N);
    end
else
    t_init = linspace(0,T0,N);
    % traj_init.x = PPTrajectory(foh([0 T0/2 T0],[x0, xm, x1]));
    traj_init.x = PPTrajectory(foh([0 T0],[x0, x1]));
    traj_init.u = PPTrajectory(zoh(t_init,0.1*randn(nu,N)));
end

T_span = [.75 T0];

optimoptions.joint_limit_collisions = false;
optimoptions.s_weight = 100;
traj_opt = VariationalTrajectoryOptimization(p,N,T_span,optimoptions);
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q0),1);  
% traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(qm),7);
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q1),N);
traj_opt = traj_opt.addVelocityConstraint(ConstantConstraint(zeros(nv,1)),1);

[q_lb, q_ub] = getJointLimits(p);
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(q_lb,q_ub),2:N-1);
z_ub = q0(3) + 0.05;
z_lb = q0(3) - 0.05;
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(z_lb,z_ub),2:N-1, 3);


state_cost = Point(getStateFrame(p),ones(nx,1));
state_cost.base_x = 0;
state_cost.base_y = 0;
state_cost.base_z = 0;
state_cost.base_pitch = 1;
state_cost.base_roll = 5;
state_cost.base_yaw = 5;
state_cost.front_left_hip_roll = 5;
state_cost.front_left_hip_pitch = 1;
state_cost.front_right_hip_roll = 5;
state_cost.front_right_hip_pitch = 5;
state_cost.back_left_hip_roll = 5;
state_cost.back_left_hip_pitch = 1;
state_cost.back_right_hip_roll = 5;
state_cost.back_right_hip_pitch = 1;
state_cost = double(state_cost);
Q = diag(state_cost); 
%Q = zeros(nx);

% traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(qf(1:6)),N);
%traj_opt = traj_opt.addInputConstraint(ConstantConstraint(zeros(3,1)),1:N-1);
%traj_opt = traj_opt.addStateConstraint(ConstantConstraint(qm(1:6)),8);
% traj_opt = traj_opt.addStateConstraint(ConstantConstraint(qf(1:6)),N);
% traj_opt = traj_opt.addPositionConstraint(periodic_constraint,{[1 N]});

traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',10000);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',1000000);
traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',1000);

traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);

tic
[xtraj,utraj,ctraj,btraj,psitraj,etatraj,~,~,straj,z,F,info] = traj_opt.solveTraj(t_init,traj_init);
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
        pause(h(1)/2);
    end
end

end
