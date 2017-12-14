function [hamr,xtraj,utraj,ctraj,btraj,...
    psitraj,etatraj,jltraj, kltraj, straj, ...
    z,F,info,infeasible_constraint_name] = HAMRVariationalPeriodicTrajOptWS(traj_init)

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
optimoptions.s_weight = 500;
optimoptions.joint_limit_collisions = false;
% optimoptions.

% ---- Initial Guess ----%
tt = traj_init.t;
T = tt(end);
N = 31;
t_init = linspace(0, T, N);
T_span = [T T];
xx0 = traj_init.x.eval(t_init(1));
xxf = traj_init.x.eval(t_init(N));

traj_opt = VariationalTrajectoryOptimization(hamr,N,T_span,optimoptions);

% -- Costs ---%
traj_opt = traj_opt.addFinalCost(@final_cost_fun);
% traj_opt = traj_opt.addRunningCost(@running_cost_fun);
% traj_opt = traj_opt.addRunningCost(@lift_cost_fun);
traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);

% -----Periodic Constraint-----

% Position
Qperiodic = [-eye(nq-1), eye(nq-1)];
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

% Add Initial/Final State Constraints
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(xx0(1:2)),1, 1:2);
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(xxf(1)-1, Inf),N,1);

% Add joint limit constraints
jlmin = hamr.joint_limit_min; jlmax = hamr.joint_limit_max;
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(jlmin,jlmax),1:N);

% Add Periodic Constraints
traj_opt = traj_opt.addPositionConstraint(periodic_cnst_pos,{[1 N]}, 2:nq);
traj_opt = traj_opt.addConstraint(periodic_cnst_vel, periodic_inds);

% Input constraints
traj_opt = traj_opt.addInputConstraint(BoundingBoxConstraint(umin, umax),1:N-1);

% ----- Solver options -----
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

%     function [f, df] = tracking_cost(x, xd, phid)
%         dim = 3;            % RigidBodyPose
%
%         Q1 = eye(dim);
%         Q2 = eye(numel(phid));
%
%         kinsol = doKinematics(hamr, x, 0*x);
%         [phi,~,~,~,~,~,~,~,n] = hamr.contactConstraints(kinsol);
%
%         f = (1/2)*(x(1:dim)-xd(1:dim))'*Q1*(x(1:dim)-xd(1:dim)) + (1/2)*(phi-phid)'*Q2*(phi-phid);
%         df = [(x(1:dim)-xd(1:dim))'*Q1, zeros(1, nq-dim)] + (phi-phid)'*Q2*n;
%         fprintf('Objective: %f \r', f)
%     end
%



    function [f,df] = final_cost_fun(tf,x)
        a = 1;
        Q = a*eye(nq+nv-1);
        f = -a*x(1) + (1/2)*(x(2:end)-xxf(2:end))'*Q'*(x(2:end)-xxf(2:end));
        df = [0, -a, (x(2:end)-xxf(2:end))'*Q];
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

    function [f,df] = running_cost_fun(h,x,u)
        R = (2/(FlimS+FLimL))^2*eye(nu);
        g = (1/2)*u'*R*u;
        f = h*g;
        df = [g, zeros(1,nx), h*u'*R];
    end


    function [f, df] = lift_cost_fun(~, x, ~)
        q = x(1:nq);
        qd = x(nq+(1:nv));
        
        kinsol = hamr.doKinematics(q, qd);
        [phi,~,~,~,~,~,~,~,n] = hamr.contactConstraints(kinsol);
        
        a = -0.2;
        lc = a*ones(1,numel(phi));
        
        f = lc*phi;
        df = [0, lc*n, zeros(1,nv), zeros(1, nu)];
        
    end
end
