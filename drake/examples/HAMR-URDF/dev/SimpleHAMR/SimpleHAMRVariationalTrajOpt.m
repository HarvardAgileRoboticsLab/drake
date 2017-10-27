function [hamr,xtraj,utraj,ctraj,btraj,...
    psitraj,etatraj,jltraj, kltraj, straj, ...
    z,F,info,infeasible_constraint_name] = SimpleHAMRVariationalTrajOpt()

% file
urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'dev', 'SimpleHAMR', ...
    'urdf', 'HAMRSimple_scaled.urdf');

% options
options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = true;
options.collision = true;

hamr = HAMRSimpleRBM(urdf,options);
v = hamr.constructVisualizer();

% state/input dimenisons
nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nx = nq+nv;
nu = hamr.getNumInputs();

% --- Set Input limits ---
ulim = 0.15;                 % set max torque
umin = -ulim*ones(nu,1);
umax = ulim*ones(nu, 1);

% --- Initialize TrajOpt---
optimoptions.s_weight = 200;
optimoptions.joint_limit_collisions = false; 
optimoptions.add_ccost = false; 

% ---- Initial Guess ----
traj_init = load('TrajOpt_26-Oct-2017_5');
t_init = traj_init.xtraj.getBreaks(); 
T = t_init(end); 
N = numel(t_init); 
x0 = traj_init.xtraj.eval(t_init(1));
q0 = x0(1:nq); v0 = 0*q0;
x1 = traj_init.xtraj.eval(T);
q1 = x1(1:nq); 

T_span = [T T];
traj_opt = VariationalTrajectoryOptimization(hamr,N,T_span,optimoptions);

init_traj.x = traj_init.xtraj;
init_traj.u = traj_init.utraj;
init_traj.c = traj_init.ctraj;
init_traj.b = traj_init.btraj;
init_traj.eta = traj_init.etatraj;
init_traj.psi = traj_init.psitraj;
init_traj.s = traj_init.straj;


% T = 200;
% N = 21;
% x0 = hamr.getInitialState();
% x1 = x0; x1(1) = x1(1) +20; q1 = x1(1:nq); 
% t_init = linspace(0,T,N);

% traj_init.x = PPTrajectory(foh([0 T],[x0, x1]));
% traj_init.u = PPTrajectory(zoh(t_init,0.1*randn(nu,N)));
% traj_init.c = PPTrajectory(zoh(t_init,0.001*randn(traj_opt.nC,N)));
% traj_init.b = PPTrajectory(zoh(t_init,0.001*randn(traj_opt.nC*traj_opt.nD,N)));
% traj_init.psi = PPTrajectory(zoh(t_init,0.001*randn(traj_opt.nC,N)));
% traj_init.eta =  PPTrajectory(zoh(t_init,0.001*randn(traj_opt.nC*traj_opt.nD,N)));

% -- State Costs ---%
state_cost = Point(getStateFrame(hamr),ones(nx,1));

state_cost.base_x = 0;
state_cost.base_y = 0;
state_cost.base_z = 1;                      % allow +/- deviations for z 

% allow +/- pi/4 rpy
state_cost.base_pitch = (4/pi)^2;           
state_cost.base_roll = (4/pi)^2;            
state_cost.base_yaw = (4/pi)^2;            

% joint limits to +/- pi/9
state_cost.FL_lift = 0; %(9/pi)^2;
% state_cost.FL_swing = (9/pi)^2;
state_cost.RL_lift = 0; %(9/pi)^2;
% state_cost.RL_swing = (9/pi)^2;
state_cost.FR_lift = 0; %(9/pi)^2;
% state_cost.FR_swing = (9/pi)^2;
state_cost.RR_lift = 0; %(9/pi)^2;
% state_cost.RR_swing = (9/pi)^2;

state_cost = double(state_cost);
state_cost(nq+(1:nv)) = 0; 
Q = diag(state_cost); 

% --- Cost Functions ---%
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
traj_opt = traj_opt.addRunningCost(@lift_cost_fun);


% -- State Constraints ---%

[qmin, qmax] = hamr.getJointLimits(); 
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(qmin, qmax),1:N);

zlb = q0(3) -1;
zub = q0(3) + 1;
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(zlb, zub),1:N, 3);

traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q0),1);
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q1(1:2)), N, 1:2);
traj_opt = traj_opt.addVelocityConstraint(ConstantConstraint(v0),1);

% Input Constraint 
traj_opt = traj_opt.addInputConstraint(BoundingBoxConstraint(umin, umax),1:N-1);

traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);


% Solver options
traj_opt = traj_opt.setSolver('snopt');
traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',10000);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',5000000);
traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',1000);

% traj_opt = traj_opt.setSolverOptions('snopt','MajorOptimalityTolerance',1e-4);
% traj_opt = traj_opt.setSolverOptions('snopt','MinorOptimalityTolerance',1e-4);
% traj_opt = traj_opt.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-4);
% traj_opt = traj_opt.setSolverOptions('snopt','MinorFeasibilityTolerance',1e-4);
% traj_opt = traj_opt.setSolverOptions('snopt','constraint_err_tol',1e-4);

disp('Solving...')
tic
[xtraj,utraj,ctraj,btraj,psitraj,etatraj,jltraj,kltraj,straj ...
    ,z,F,info,infeasible_constraint_name] = solveTraj(traj_opt,t_init,init_traj);
toc
%
    function [f,df] = running_cost_fun(h,x,u)
        R = (1/ulim)^2*eye(nu);
        g = (1/2)*(x-x1)'*Q*(x-x1) + (1/2)*u'*R*u;
        f = h*g;
        df = [g, h*(x-x1)'*Q, h*u'*R];
    end
   

    function [f, df] = lift_cost_fun(h, x, u) 
        q = x(1:nq);
        qd = x(nq+(1:nv)); 
%         lift_cost = Point(getStateFrame(hamr),zeros(nx,1));
%         lift_cost.FL_lift = -0.1; 
%         lift_cost.RL_lift = -0.1; 
%         lift_cost.FR_lift = 0.1; 
%         lift_cost.RR_lift = 0.1; 
%         lift_cost = double(lift_cost); 
        kinsol = hamr.doKinematics(q, qd, struct('compute_gradients',true)); 
        [phi,~,~,~,~,~,~,~,n,~,~,~] = hamr.contactConstraints(kinsol);
        lift_cost = 0.05*ones(size(phi)); 
        f = -lift_cost'*phi; 
        df = [0, -lift_cost'*n, zeros(1,nv), zeros(1, nu)]; 

    end

    function [f,df] = final_cost_fun(tf,x)
        a = 0.1;
        f = -a*x(1);
        df = zeros(1, nx+1);
        df(2) = -a;
    end

    function displayTraj(h,x,u)
        disp('Displaying Trajectory...')
        h = h/1e3;
        ts = [0;cumsum(h)];
        for i=1:length(ts)
            v.drawWrapper(0,x(:,i));
            pause(h(1));
        end
        
    end
end
