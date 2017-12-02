function [hamr,xtraj,utraj,ctraj,btraj,...
    psitraj,etatraj,jltraj, kltraj, straj, ...
    z,F,info,infeasible_constraint_name] = SimpleHAMRVariationalTrajOpt(save_dir)

% file
urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'urdf', 'HAMRSimple_scaled.urdf');

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
Vlim = 100;                     % set max torque
V2Tau = 0.0013; 
umin = -Vlim*V2Tau*ones(nu,1);
umax = Vlim*V2Tau*ones(nu, 1);

% --- Initialize TrajOpt---
optimoptions.s_weight = 1200;
optimoptions.joint_limit_collisions = false; 
optimoptions.add_ccost = true; 

% ---- Initial Guess ----
fname = 'TrajOpt_MovingBody_SimpleSprings7';
traj0 = load([save_dir, fname]);
t_init = traj0.xtraj.getBreaks(); 
N = 1*(numel(t_init)-1)+1; 
if N < numel(t_init)
    T = t_init(N); 
else 
    T = t_init(end);
end

t_init = linspace(0, T, N); 
x0 = traj0.xtraj.eval(t_init(1));
x1 = traj0.xtraj.eval(t_init(N)); 
q1 = x1(1:nq); 
traj_init.x = traj0.xtraj;
traj_init.u = traj0.utraj;
traj_init.c = traj0.ctraj;
traj_init.b = traj0.btraj;
traj_init.eta = traj0.etatraj;
traj_init.psi = traj0.psitraj;
% traj_init.s = traj0.straj;

% -- Initialize Traj Opt ---% 

% T = 100;
% N = 21;
% x0 = hamr.getInitialState();
% x1 = x0; x1(1) = x1(1); q1 = x1(1:nq);

T_span = [T T];
traj_opt = VariationalTrajectoryOptimization(hamr,N,T_span,optimoptions);
 
% t_init = linspace(0,T,N);
% traj_init.x = PPTrajectory(foh([0 T],[x0, x1]));
% traj_init.u = PPTrajectory(zoh(t_init, 0.001*randn(nu,N)));
% traj_init.c = PPTrajectory(zoh(t_init,0.001*randn(traj_opt.nC,N)));
% traj_init.b = PPTrajectory(zoh(t_init,0.001*randn(traj_opt.nC*traj_opt.nD,N)));
% traj_init.psi = PPTrajectory(zoh(t_init,0.001*randn(traj_opt.nC,N)));
% traj_init.eta =  PPTrajectory(zoh(t_init,0.001*randn(traj_opt.nC*traj_opt.nD,N)));


% -- State Costs ---%
state_cost = Point(getStateFrame(hamr),ones(nx,1));

state_cost.base_x = 0;
state_cost.base_y = 0;
state_cost.base_z = 1;                       

% allow +/- pi/4 rpy
state_cost.base_pitch = (8/pi)^2;           
state_cost.base_roll = (8/pi)^2;            
state_cost.base_yaw = (8/pi)^2;            

% joint limits to +/- pi/9
state_cost.FL_lift = 0; %(9/pi)^2;
state_cost.FL_swing = 0; %(9/pi)^2;
state_cost.RL_lift = 0; %(9/pi)^2;
state_cost.RL_swing = 0; %(9/pi)^2;
state_cost.FR_lift = 0; %(9/pi)^2;
state_cost.FR_swing = 0; %(9/pi)^2;
state_cost.RR_lift = 0; %(9/pi)^2;
state_cost.RR_swing = 0; %(9/pi)^2;

state_cost = double(state_cost);
state_cost(nq+(1:nv)) = 0; 
Q = diag(state_cost); 

% --- Cost Functions ---%
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
traj_opt = traj_opt.addRunningCost(@lift_cost_fun);
traj_opt = traj_opt.addFinalCost(@final_cost_fun); 

% -- State Constraints ---%
[qmin, qmax] = hamr.getJointLimits(); 
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(qmin, qmax),1:N);

zlb = x0(3) - 1;
zub = x0(3) + 1;
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(zlb, zub),1:N, 3);
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(x0(1:nq)),1);
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q1(1:2)), N, 1:2);
traj_opt = traj_opt.addVelocityConstraint(ConstantConstraint(x0(nq+(1:nv))),1);

% Input Constraint 
traj_opt = traj_opt.addInputConstraint(BoundingBoxConstraint(umin, umax),1:N-1);
traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);


% Solver options
traj_opt = traj_opt.setSolver('snopt');
traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',10000);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',5000000);
traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',1000);

traj_opt = traj_opt.setSolverOptions('snopt','MajorOptimalityTolerance',1e-6);
traj_opt = traj_opt.setSolverOptions('snopt','MinorOptimalityTolerance',1e-6);
traj_opt = traj_opt.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-6);
traj_opt = traj_opt.setSolverOptions('snopt','MinorFeasibilityTolerance',1e-6);
traj_opt = traj_opt.setSolverOptions('snopt','constraint_err_tol',1e-6);

disp('Solving...')
tic
[xtraj,utraj,ctraj,btraj,psitraj,etatraj,jltraj,kltraj,straj ...
    ,z,F,info,infeasible_constraint_name] = solveTraj(traj_opt,t_init,traj_init);
toc

    function [f,df] = running_cost_fun(h,x,u)
        ulim = Vlim*V2Tau; 
        R = (1/ulim)^2*eye(nu);
        g = (1/2)*(x-x1)'*Q*(x-x1) + (1/2)*u'*R*u;
        f = h*g;
        df = [g, h*(x-x1)'*Q, h*u'*R];
    end

    function [f, df] = lift_cost_fun(h, x, u)
        xin = [h; x; u]; 
        [f, df] = lift_cost(xin);        
%         
%         df_fd = zeros(size(df));
%         step = 1e-6;
%         dxin = step*eye(length(xin));
%         for k = 1:length(xin)
%             df_fd(:,k) = (lift_cost(xin+dxin(:,k)) - lift_cost(xin-dxin(:,k)))/(2*step);
%         end
%         
%         disp('Lift Cost Derivative Error:');
%         disp(max(abs(df_fd(:)-df(:))));        
    end            
   

    function [f, df] = lift_cost(xin) 
        h = xin(1); 
        q = xin(1+(1:nq));
        qd = xin(1+(nq+(1:nv))); 
        
        kinsol = hamr.doKinematics(q, qd, struct('compute_gradients', true)); 
        [phi,~,~,~,~,~,~,~,n] = hamr.contactConstraints(kinsol);        
 
%         a = -1; 
%         lc = a*ones(1,numel(phi)); 
        lc = [-1, -0.1, -1, -0.1]; 
%         lc([8, 10]) = -10; 
%         lc([12, 14]) = 10; 

        f = lc*phi; 
        df = [0, lc*n, zeros(1,nv), zeros(1, nu)]; 

    end

    function [f,df] = final_cost_fun(tf,x)
        a = 10;
        f = a*x(5)^2;
        df = zeros(1, nx+1);
        df(6) = 2*a*x(5);
    end

    function displayTraj(h,x,u)
        disp('Displaying Trajectory...')
        h = h/1e3;
        ts = [0;cumsum(h)];
        for i=1:length(ts)
            v.drawWrapper(0,x(:,i));
            pause(5*h(1));
        end
        
    end
end
