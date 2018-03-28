function [hamr,xtraj,utraj,ctraj,btraj,...
    psitraj,etatraj,jltraj, kltraj, straj, ...
    z,F,info,infeasible_constraint_name] = HAMRVariationalTrajOptWS(traj_init)

% file
urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2.urdf');

% options
options.terrain = RigidBodyFlatTerrain();
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
nact = hamr.getNumActuatedDOF(); 

% --- Set Input limits ---
FlimS = 0.30;                 % set max force
FLimL = 0.30;
umin = -reshape([FlimS; FLimL]*ones(1, nu/2), nu, 1);
umax =  reshape([FlimS; FLimL]*ones(1, nu/2), nu, 1);
hamr = hamr.setInputLimits(umin, umax);

% --- Initialize TrajOpt---
optimoptions.s_weight = 200;
optimoptions.joint_limit_collisions = false;

% ---- Initial Guess ----%
tt = traj_init.t;
T = tt(end);
N = 21;
t_init = linspace(0, T, N);
T_span = [T T];
xx0 = traj_init.x.eval(t_init(1));
xxf = traj_init.x.eval(t_init(N));

% xxf(1:2) = xxf(1:2) - xx0(1:2);
% xx0(1:2) = 0;

traj_opt = VariationalTrajectoryOptimization(hamr,N,T_span,optimoptions);

% -- Costs ---%
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
% nvarsCOTMech = (nu+nact)*N + 1; 
% COTMech = FunctionHandleObjective(nvarsCOTMech,@(x)tracking_cost(x),1), ...
%         {traj_opt.x_inds(,ind)

traj_opt = traj_opt.addRunningCost(@lift_cost_fun);

% -- Constraints ---%
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(xx0(1:nq, 1)),1);
traj_opt = traj_opt.addVelocityConstraint(ConstantConstraint(xx0(nv+(1:nq), 1)),1);
traj_opt =  traj_opt.addPositionConstraint(ConstantConstraint(xxf(1:nq, 1)),N);
traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);

jlmin = hamr.joint_limit_min; %jlmin(3) = q0(3) - 1.5;
jlmax = hamr.joint_limit_max; %jlmax(3) = q0(3) + 1.5;

traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(jlmin,jlmax),2:N);
traj_opt = traj_opt.addInputConstraint(BoundingBoxConstraint(umin, umax),1:N-1);

% Solver options
traj_opt = traj_opt.setSolver('snopt');
traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',10000);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',10000000);
traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',1000);
% traj_opt = traj_opt.setSolverOptions('snopt','print','outputlog.txt');

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

    function [f,df] = running_cost_fun(h,x,u)
        R = (2/(FlimS+FLimL))^2*eye(nu);
        g = (1/2)*u'*R*u;
        f = h*g;
        df = [g, zeros(1,nx), h*u'*R];
    end


    function [f, df] = lift_cost_fun(xin)
        h = xin(1);
        q = xin(1+(1:nq));
        qd = xin(1+(nq+(1:nv)));
        
        kinsol = hamr.doKinematics(q, qd);
        [phi,~,~,~,~,~,~,~,n] = hamr.contactConstraints(kinsol);
        
        a = -0.2;
        lc = a*ones(1,numel(phi));
        
        f = lc*phi;
        df = [0, lc*n, zeros(1,nv), zeros(1, nu)];
        
    end

    function [f, df] = tracking_cost(x, xd, phid)
        dim = 3;            % RigidBodyPose
        
        Q1 = eye(dim);
        Q2 = eye(numel(phid));
        
        kinsol = doKinematics(hamr, x, 0*x);
        [phi,~,~,~,~,~,~,~,n] = hamr.contactConstraints(kinsol);
        
        f = (1/2)*(x(1:dim)-xd(1:dim))'*Q1*(x(1:dim)-xd(1:dim)) + (1/2)*(phi-phid)'*Q2*(phi-phid);
        df = [(x(1:dim)-xd(1:dim))'*Q1, zeros(1, nq-dim)] + (phi-phid)'*Q2*n;
        fprintf('Objective: %f \r', f)
    end
%
%
%     function [f,df] = final_cost_fun(tf,x)
%         a = 1;
%         f = -a*x(1);
%         df = zeros(1, nx+1);
%         df(2) = -a;
%     end

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
