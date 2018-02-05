function [hamr,xtraj,utraj,ctraj,btraj,...
    psitraj,etatraj,jltraj, kltraj, straj, ...
    z,F,info,infeasible_constraint_name] = HAMRVariationalStepUp(traj_init)
%% Set up robot
urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2.urdf');

% options
options.terrain = RigidBodyFlatTerrain;
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = true;
options.collision = true;

hamr = HamrRBM(urdf,options);
v = hamr.constructVisualizer();

% state/input dimenisons
nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nx = nq+nv;
nu = hamr.getNumInputs();
% act_dof = hamr.getActuatedJoints(); 

% --- Set Input limits ---
FlimS = 0.28;                 % set max force
FLimL = 0.36;
umin = -reshape([FlimS; FLimL]*ones(1, nu/2), nu, 1);
umax =  reshape([FlimS; FLimL]*ones(1, nu/2), nu, 1);
hamr = hamr.setInputLimits(umin, umax);

%% Set up Traj Opt

% --- Initialize TrajOpt---
optimoptions.s_weight = 200;
optimoptions.joint_limit_collisions = false;

% ---- Initial Guess ----%
tt = traj_init.t;
T = tt(end);
N = 31;
t_init = linspace(0, T, N);
x0 = traj_init.x.eval(t_init(1));
x1 = traj_init.x.eval(t_init(end));
t_init = linspace(0, T, N);

T_span = [T T];

% Initialzie optimization
traj_opt = VariationalTrajectoryOptimization3(hamr,N,T_span,optimoptions);

% -- Costs ---%
% traj_opt = traj_opt.addCost(

% try to get smooth u
% for ind = 1:N-2
%     traj_opt = traj_opt.addCost(FunctionHandleObjective(2*nu,@(u1, u2)udiff_cost(u1, u2)), ...
%         {traj_opt.u_inds(:,ind); traj_opt.u_inds(:,ind+1)});
% end

traj_opt = traj_opt.addFinalCost(@final_cost); 
traj_opt = traj_opt.addRunningCost(@running_cost);
traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);

% -- Constraints ---%

% joint limit constraints
[qmin, qmax] = getJointLimits(hamr);
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(qmin,qmax),1:N);

% x-constriant
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(x0(1:nq)),1);
traj_opt = traj_opt.addVelocityConstraint(ConstantConstraint(x0(nq+(1:nq))),1);

% Input constraints
traj_opt = traj_opt.addInputConstraint(BoundingBoxConstraint(umin, umax),1:N-1);

% StateCosts
state_cost = Point(getStateFrame(hamr),ones(nx,1));
state_cost.base_x = 1; 
state_cost.base_y = 1;
state_cost.base_z = 10;
state_cost.base_pitch = (10/pi)^2;
state_cost.base_roll = (10/pi)^2;
state_cost.base_yaw = (10/pi)^2;
state_cost.base_zdot = 10; 
state_cost = double(state_cost);

% state_cost(7:end) = 0;
% state_cost(act_dof(1:2:end)) = 5; 
Q = diag(state_cost);



% ------- Solver options ----------
traj_opt = traj_opt.setSolver('snopt');
traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',10000);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',10000000);
traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',1000);

traj_opt = traj_opt.setSolverOptions('snopt','MajorOptimalityTolerance',1e-5);
traj_opt = traj_opt.setSolverOptions('snopt','MinorOptimalityTolerance',1e-5);
traj_opt = traj_opt.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-5);
traj_opt = traj_opt.setSolverOptions('snopt','MinorFeasibilityTolerance',1e-5);
traj_opt = traj_opt.setSolverOptions('snopt','constraint_err_tol',1e-5);


disp('Solving...')
tic
[xtraj,utraj,ctraj,btraj, psitraj,etatraj,jltraj, kltraj,straj, ...
    z,F,info,infeasible_constraint_name] = traj_opt.solveTraj(t_init,traj_init);
toc
% 
    function [c,dc] = final_cost(tf,xf)
        c = (1/2)*(xf-x1)'*Q*(xf-x1);
        dc = [0, (xf-x1)'*Q];
        fprintf('Final Cost %f \r', c);
    end

    function [f,df] = running_cost(h, x, u)     
        r = ones(1, numel(u)); 
        r(2:2:end) = 0; 
        R = diag(r);         
        g =(1/2)*u'*R*u;
        f = h*g;
        df = [g, zeros(1, numel(x)), h*u'*R];
    end

% 
    function [f, df] = udiff_cost(u1, u2)
        udiff = u1-u2;
        f = 0.5*(udiff'*udiff);
        I = eye(length(u1));
        df = [udiff'*I,-udiff'*I];
    end

    function displayTraj(h,x,u)
        disp('Displaying Trajectory...')
        h = h/1e3;
        ts = [0;cumsum(h)];
        for i=1:length(ts)
            v.drawWrapper(0,x(:,i));
            pause(10*h(1));
        end
        
    end
end
