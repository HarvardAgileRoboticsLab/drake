function [xtraj,utraj,ctraj,btraj,psitraj,etatraj,straj] = runVariational(x_guess,u_guess,c_guess,b_guess,psi_guess,eta_guess,s_guess)

options=struct();
options.terrain = RigidBodyFlatTerrain();
options.use_bullet = false;
options.ignore_self_collisions = true;
options.multiple_contacts = false;
options.active_collision_options.terrain_only = false;
options.floor_off = false;
r = KukaArm(options);

nq = r.getNumPositions();
nv = r.getNumVelocities();
nx = nq+nv;
nu = r.getNumInputs();

v=r.constructVisualizer;

q0 = [-1.57;-1.1;0;1.57;0.0;0.1;0;0.08; ...
      0;0.66;0.03;0;0;0];
x0 = [q0;zeros(nq,1)];
v.draw(0,x0);

q1 = q0;
q1(11) = q0(11)+0.02;
x1 = [q1;zeros(nq,1)];

u0 = r.findTrim(q0);
u0(8) = -10;

T0 = 1;
N = 5;
T_span = [.5 T0];

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
    if length(x_guess.getBreaks()) == N
        t_init = x_guess.getBreaks();
    else
        t_init = linspace(x_guess.tspan(1),x_guess.tspan(end),N);
    end
else
    t_init = linspace(0,T0,N);
    traj_init.x = PPTrajectory(foh([0 T0],[x0, x0]));
    traj_init.u = PPTrajectory(zoh([0 T0],[u0, u0]));
    traj_init.c = PPTrajectory(zoh([0 T0],[[10*9.81; 0; 0; 0], [10*9.81; 0; 0; 0]]));
end

options.lambda_weight = 0.01;
options.s_weight = 1000;
options.s0 = .01;
options.s_max = .01;
options.s_min = 1e-6;

traj_opt = VariationalTrajectoryOptimization(r,N,T_span,options);
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q0),1); 
%traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q0),N); 
%traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q1(9:11)),1:N,9:11);
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q1(9:11)),N,9:11);

[q_lb, q_ub] = getJointLimits(r);
q_lb = max([q_lb, q0-0.2*ones(14,1)]')';
q_ub = min([q_ub, q0+0.2*ones(14,1)]')';
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(q_lb,q_ub),1:N);

% ub_N = q1;
% ub_N(1:8) = q_ub(1:8);
% lb_N = q1;
% lb_N(1:8) = q_lb(1:8);
% traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(lb_N,ub_N),N);
% traj_opt = traj_opt.addVelocityConstraint(ConstantConstraint(zeros(nv,1)),1);

% state_cost = Point(getStateFrame(r),ones(nx,1));
% state_cost.base_x = 10;
% state_cost.base_y = 10;
% state_cost.base_z = 10;
% state_cost = double(state_cost);
Q = blkdiag(1*eye(8),0*eye(6),10*eye(14)); %diag(state_cost); 
R = 1*eye(nu);

%traj_opt = traj_opt.setSolver('fmincon');
% traj_opt = traj_opt.setSolver('ipopt');

traj_opt = traj_opt.setSolverOptions('snopt','majorfeasibilitytolerance',1e-4);
traj_opt = traj_opt.setSolverOptions('snopt','minorfeasibilitytolerance',1e-4);
traj_opt = traj_opt.setSolverOptions('snopt','minoroptimalitytolerance',1e-4);
traj_opt = traj_opt.setSolverOptions('snopt','majoroptimalitytolerance',1e-4);
traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',100000);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',2000000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',10000000);
traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',10000);

traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);

tic
[xtraj,utraj,ctraj,btraj,psitraj,etatraj,~,~,straj,z,F,info] = traj_opt.solveTraj(t_init,traj_init);
toc

v.playback(xtraj,struct('slider',true));

function [f,df] = running_cost_fun(h,x,u)
  g = (1/2)*(x-x0)'*Q*(x-x0) + (1/2)*(u-u0)'*R*(u-u0);
  f = h*g;
  df = [g, h*(x-x0)'*Q, h*(u-u0)'*R];
end

function displayTraj(h,x,u)
  
  ts = [0;cumsum(h)];
  for i=1:length(ts)
    v.drawWrapper(0,x(:,i));
    pause(h(1)/4);
  end
   
end

end