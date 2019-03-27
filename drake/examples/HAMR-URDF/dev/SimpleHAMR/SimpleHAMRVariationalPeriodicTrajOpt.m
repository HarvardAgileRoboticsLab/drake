function [hamr,x,u,c,b, psi,eta,jl, kl, s, z,F,info, ...
    infeasible_constraint_name] = SimpleHAMRVariationalPeriodicTrajOpt(load_path, params)

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
ulim = 0.13;                 % set max force
umin = -ulim*ones(nu,1);
umax = ulim*ones(nu, 1);

% --- Initialize TrajOpt---
optimoptions.s_weight = params.s_weight;
optimoptions.add_ccost = true;

% parameters
T = 10;  % ms
DX = 10; % mm
xi = hamr.getInitialState(); xi(nq+1) = DX/T;
xf = xi; xf(1) = DX; xf(nq+1) = DX/T;
N = 31;

if isempty(load_path)
    t_init = linspace(0, T, N);
    traj_init.x = PPTrajectory(foh([0 T],[xi, xf]));
    traj_init.u = PPTrajectory(zoh(t_init,0.1*randn(nu,N)));
else
    traj_init = load(load_path);
    t_init  = traj_init.x.getBreaks();
    traj_init = rmfield(traj_init, {'jl', 'kl'});
    traj_init.s = 0*traj_init.s + 1; 
    
end

T_span = [0.5*T 3*T];
traj_opt = VariationalTrajectoryOptimization(hamr,N,T_span,optimoptions);

% -- Costs ---%
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
traj_opt = traj_opt.addRunningCost(@foot_height_fun); 
traj_opt = traj_opt.addFinalCost(@final_cost_fun);

% traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);

% -----Periodic Constraint-----%

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
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(xi(1:6)),1,1:6);
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(xf(1:6)-1, xf(1:6)+1),N,1:6);

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

traj_opt = traj_opt.setSolverOptions('snopt','MajorOptimalityTolerance',1e-5);
traj_opt = traj_opt.setSolverOptions('snopt','MinorOptimalityTolerance',1e-5);
traj_opt = traj_opt.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-5);
traj_opt = traj_opt.setSolverOptions('snopt','MinorFeasibilityTolerance',1e-5);
traj_opt = traj_opt.setSolverOptions('snopt','constraint_err_tol',1e-5);


disp('Solving...')
tic
[x,u,c,b,psi,eta,jl,kl,s,z,F,info,infeasible_constraint_name] ...
    = solveTraj(traj_opt,t_init,traj_init);
toc

    function [f,df] = running_cost_fun(h,x,u)
        R = (1/ulim)^2*eye(nu);
        g = (1/2)*u'*R*u;
        f = h*g;
        df = [g, zeros(1,nx), h*u'*R];
    end

    function [f,df] = final_cost_fun(tf,x)
        Q = diag([1, 1, 1, 10, 10, 10]);
        f = 0.5*(x(1:6) - xf(1:6))'*Q*(x(1:6) - xf(1:6));
        df = [0, (x(1:6) - xf(1:6))'*Q, zeros(1, (nq+nv)-6)];
    end

    function [f,df] = foot_height_fun(h,x,u)
        q = x(1:nq);
        
        [phi,~,~,~,~,~,~,~,n] = hamr.contactConstraints(q,false,struct('terrain_only',true));
        phi0 = [.25;.25;.25;.25];
        K = 20;
        I = find(phi < phi0);
        f = K*(phi(I) - phi0(I))'*(phi(I) - phi0(I));
        % phi: 2x1
        % n: 2xnq
        df = [0 2*K*(phi(I)-phi0(I))'*n(I,:) zeros(1,nv+nu)];
        

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

    function [f,df] = periodic_constraint_fun(h0,q0,q1,u0,c0,b0,jl0,kl0, ...
            hNm1,qNm1,qN,uNm1,klNm1)
        
        xin = [h0;q0;q1;u0;c0;b0;jl0;kl0; ...
            hNm1;qNm1;qN;uNm1;klNm1];
        [f,df] = periodic_constraint(xin);
        
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
end
