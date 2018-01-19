function [hamr,xtraj,utraj,ctraj,btraj,...
    psitraj,etatraj,jltraj, kltraj, straj, ...
    z,F,info,infeasible_constraint_name] = HAMRVariationalStepUp(varargin)
%% Set up Steps
% w  = 100;
% l = 50;
% h = 5;

% boxes = [0, 0, l, w, h];


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
act_dof = hamr.getActuatedJoints(); 

% --- Set Input limits ---
FlimS = 0.3;                 % set max force
FLimL = 0.3;
umin = -reshape([FlimS; FLimL]*ones(1, nu/2), nu, 1);
umax =  reshape([FlimS; FLimL]*ones(1, nu/2), nu, 1);
hamr = hamr.setInputLimits(umin, umax);

%% Set up Traj Opt

% --- Initialize TrajOpt---
optimoptions.s_weight = 600;
optimoptions.joint_limit_collisions = false;

% ---- Initial Guess ----%

if ~isempty(varargin)
    traj0 = varargin{1};
    t_init = traj0.xtraj.getBreaks();
    N = numel(t_init);
    T = t_init(end);
    x0 = traj0.xtraj.eval(t_init(1));
    x1 = traj0.xtraj.eval(t_init(end));
else
    T = 100;
    N = 41;
    t_init = linspace(0, T, N);
    x0 = hamr.getInitialState(); x0(3) = 12.59;
    x1 = x0; 
    
end


T_span = [0.8*T 1.2*T];

% Initialzie optimization
traj_opt = VariationalTrajectoryOptimization(hamr,N,T_span,optimoptions);

if ~isempty(varargin)
    traj_init.x = traj0.xtraj; %.eval(;
    traj_init.u = traj0.utraj;
    traj_init.c = traj0.ctraj;
    traj_init.b = traj0.btraj;
    traj_init.psi = traj0.psitraj;
    traj_init.eta = traj0.etatraj;
    traj_init.kl = traj0.kltraj;
else
    traj_init.x = PPTrajectory(foh([0, T], [x0, x1]));
    traj_init.u = PPTrajectory(zoh(t_init, 0.2*randn(nu,N)));
    traj_init.c = PPTrajectory(zoh(t_init,0.001*randn(traj_opt.nC,N)));
    traj_init.b = PPTrajectory(zoh(t_init,0.001*randn(traj_opt.nC*traj_opt.nD,N)));
    traj_init.psi = PPTrajectory(zoh(t_init,0.001*randn(traj_opt.nC,N)));
    traj_init.eta =  PPTrajectory(zoh(t_init,0.001*randn(traj_opt.nC*traj_opt.nD,N)));
    traj_init.kl =  PPTrajectory(zoh(t_init,0.001*randn(traj_opt.nKL,N)));
end
% -- Costs ---%

% % try to get smooth accelerations
% for ind = 2:N-1
%     traj_opt = traj_opt.addCost(FunctionHandleObjective(2*nv,@vdiff_cost), ...
%         {traj_opt.h_inds(ind-1); traj_opt.h_inds(ind); traj_opt.x_inds(:,ind-1); ...
%         traj_opt.x_inds(:,ind); traj_opt.x_inds(:,ind+1)});
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
state_cost.base_x = 0;
state_cost.base_y = 0;
state_cost.base_z = 0;
state_cost.base_pitch = (10/pi)^2;
state_cost.base_roll = (10/pi)^2;
state_cost.base_yaw = (10/pi)^2;
state_cost = double(state_cost);
state_cost(7:end) = 0;
state_cost(act_dof(1:2:end)) = 5; 
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
        a = 20;                
        c = -a*(xf(3) - x0(3)) + (xf-x1)'*Q*(xf-x1);
        dc_dx = (xf-x1)'*Q; 
        dc_dx(3) = dc_dx(3) - a; 
        dc = [0, dc_dx];
        fprintf('Final Cost %f \r', c); 
    end

    function [f,df] = running_cost(h, x, u)        
        
        r = (1/FLimL)^2*ones(nu, 1);         
        r(2:2:end) = 0; 
        R = diag(r); 
        
        
        g =(1/2)*((x-x1)'*Q*(x-x1) + u'*R*u);
        f = h*g;
        df = [g, (x-x1)'*Q, h*u'*R];
    end

    function [c, dc] = vdiff_cost(h1, h2, q1, q2, q3)        
       
        %Take care of angle wrap-around
        vm1 = traj_opt.qdiff(q1,q2,h1);
        vm2 = traj_opt.qdiff(q2,q3,h2);
        
        a = 1; 
        vdiff = (vm1 - vm2);
        c = 0.5*a*(vdiff'*vdiff);
        I = eye(nv);
        dc = a*vdiff'*[-vm1/h1, vm2/h2, -I/h1, I/h1+I/h2, -I/h2];
        fprintf('Vdiff Cost: %f \r', c)
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
