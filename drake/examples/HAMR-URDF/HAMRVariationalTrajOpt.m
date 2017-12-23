function [hamr,xtraj,utraj,ctraj,btraj,...
    psitraj,etatraj,jltraj, kltraj, straj, ...
    z,F,info,infeasible_constraint_name] = HAMRVariationalTrajOpt()

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
T = 0.01;
T_span = [0.8*T 1.2*T];
N = 41;
t_init = linspace(0, T, N);
x0 = hamr.getInitialState(); x0(3) = 12.59; x0(1+nq) = 20/T; 
x1 = x0; x1(1) = x1(1)+20; q1 = x1(1:nq);


% Initialzie optimization
traj_opt = VariationalTrajectoryOptimization(hamr,N,T_span,optimoptions);

traj_init.x = PPTrajectory(foh([0, T], [x0, x1]));
traj_init.u = PPTrajectory(zoh(t_init, 0.1*randn(nu,N)));
traj_init.c = PPTrajectory(zoh(t_init,0.001*randn(traj_opt.nC,N)));
traj_init.b = PPTrajectory(zoh(t_init,0.001*randn(traj_opt.nC*traj_opt.nD,N)));
traj_init.psi = PPTrajectory(zoh(t_init,0.001*randn(traj_opt.nC,N)));
traj_init.eta =  PPTrajectory(zoh(t_init,0.001*randn(traj_opt.nC*traj_opt.nD,N)));
traj_init.kl =  PPTrajectory(zoh(t_init,0.001*randn(traj_opt.nKL,N)));
traj_init.s = PPTrajectory(zoh([0, T], [1, 1]));


% -- Costs ---%
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);

% -- Constraints ---%

% Position
Qperiodic = [-eye(nq-1), eye(nq-1)];        %FOR FIXED BODY
periodic_cnst_pos = LinearConstraint(zeros(nq-1,1),zeros(nq-1,1),Qperiodic);
periodic_cnst_pos = periodic_cnst_pos.setName('periodicity_pos');

% Velocity
cnstr_opts.grad_level = 1;
cnstr_opts.grad_method = 'user';
periodic_vars = 2 + 4*nq+2*nu+traj_opt.nC+traj_opt.nD*traj_opt.nC+traj_opt.nJL+2*traj_opt.nKL;
periodic_cnst_vel = FunctionHandleConstraint(zeros(nq,1), zeros(nq,1), periodic_vars, ...
    @periodic_constraint_fun, cnstr_opts);

periodic_cnst_vel = periodic_cnst_vel.setName('periodicity_vel');

periodic_inds = {traj_opt.h_inds(1); traj_opt.x_inds(:,1); traj_opt.x_inds(:,2); ...
    traj_opt.u_inds(:,1); traj_opt.c_inds(:,1); traj_opt.b_inds(:,1); traj_opt.jl_inds(:,1);
    traj_opt.kl_inds(:,1); traj_opt.h_inds(N-1); traj_opt.x_inds(:,N-1); traj_opt.x_inds(:,N); ...
    traj_opt.u_inds(:,N-1); traj_opt.kl_inds(:,N-1)};

% Add joint limit constraints
jlmin = hamr.joint_limit_min; jlmax = hamr.joint_limit_max;
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(jlmin,jlmax),1:N);

% Add x-constriant
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(x0(1)),1, 1);
% traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(x1(1), Inf) ,N, 1);

% Add Periodic Constraints
traj_opt = traj_opt.addPositionConstraint(periodic_cnst_pos,{[1 N]}, 2:nq);
traj_opt = traj_opt.addConstraint(periodic_cnst_vel, periodic_inds);

% Input constraints
traj_opt = traj_opt.addInputConstraint(BoundingBoxConstraint(umin, umax),1:N-1);


% ------- Solver options ----------
traj_opt = traj_opt.setSolver('snopt');
traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',10000);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',10000000);
traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',1000);

traj_opt = traj_opt.setSolverOptions('snopt','MajorOptimalityTolerance',1e-6);
traj_opt = traj_opt.setSolverOptions('snopt','MinorOptimalityTolerance',1e-6);
traj_opt = traj_opt.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-6);
traj_opt = traj_opt.setSolverOptions('snopt','MinorFeasibilityTolerance',1e-6);
traj_opt = traj_opt.setSolverOptions('snopt','constraint_err_tol',1e-6);


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

    function [f,df] = periodic_constraint_fun(h0,q0,q1,u0,c0,b0,jl0,kl0, ...
            hNm1,qNm1,qN,uNm1,klNm1)
        
        xin = [h0;q0;q1;u0;c0;b0;jl0;kl0; ...
            hNm1;qNm1;qN;uNm1;klNm1];
        [f,df] = periodic_constraint(xin);
        %         f
        
        %         df_fd = zeros(size(df));
        %         step = 1e-6; % sqrt(eps(max(xin)));
        %         dxin = step*eye(length(xin));
        %         for k = 1:length(xin)
        %             xin + dxin(:,k);
        %             df_fd(:,k) = (periodic_constraint(xin+dxin(:,k)) - ...
        %                 periodic_constraint(xin-dxin(:,k)))/(2*step);
        %         end
        %
        %         disp('Periodic constraint derivative error:');
        %         disp(max(abs(df_fd(:)-df(:))));
    end


    function [f, df] = periodic_constraint(xin)
        
        nC = traj_opt.nC;
        nD = traj_opt.nD;
        nJL = traj_opt.nJL;
        nKL = traj_opt.nKL;
        nQ = traj_opt.plant.getNumPositions();
        nU = traj_opt.plant.getNumInputs();
        
        % left side
        h0 = xin(1);
        q0 = xin(1+(1:nQ));
        q1 = xin(1+nQ+(1:nQ));
        u0 = xin(1+2*nQ+(1:nU));
        c0 = xin(1+2*nQ+nU+(1:nC));
        b0 = xin(1+2*nQ+nU+nC+(1:nD*nC));
        jl0 = xin(1+2*nQ+nU+nC+nD*nC+(1:nJL));
        kl0 = xin(1+2*nQ+nU+nC+nD*nC+nJL+(1:nKL));
        
        % right side
        hNm1 = xin(1+2*nQ+nU+nC+nD*nC+nJL+nKL+1);
        qNm1 = xin(1+2*nQ+nU+nC+nD*nC+nJL+nKL+1+(1:nQ));
        qN = xin(1+2*nQ+nU+nC+nD*nC+nJL+nKL+1+nQ+(1:nQ));
        uNm1 = xin(1+2*nQ+nU+nC+nD*nC+nJL+nKL+1+2*nQ+(1:nU));
        klNm1 = xin(1+2*nQ+nU+nC+nD*nC+nJL+nKL+1+2*nQ+nU+(1:nKL));
        
        [p0, dp0] = left_legendre_transform_fun(traj_opt,h0,q0,q1,u0,c0,b0,jl0,kl0);
        [pN, dpN] = right_legendre_transform_fun(traj_opt,hNm1,qNm1,qN,uNm1,klNm1);
        
        f = p0 - pN;
        
        df = zeros(nQ, numel(xin));
        df(:, 1:1+2*nQ+nU+nC+nD*nC+nJL+nKL) = dp0;
        df(:, 1+2*nQ+nU+nC+nD*nC+nJL+nKL+(1:1+2*nQ+nU+nKL)) = -dpN;
        
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
