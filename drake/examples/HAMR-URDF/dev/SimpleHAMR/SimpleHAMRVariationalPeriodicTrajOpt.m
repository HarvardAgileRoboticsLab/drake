function [hamr,xtraj,utraj,ctraj,btraj,...
    psitraj,etatraj,jltraj, kltraj, straj, ...
    z,F,info,infeasible_constraint_name] = SimpleHAMRVariationalPeriodicTrajOpt()

% file
urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'dev', 'SimpleHAMR', 'HAMRSimple_scaled.urdf');

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
ulim = 10;                 % set max force
umin = -ulim*ones(nu,1);
umax = ulim*ones(nu, 1);

% --- Initialize TrajOpt---
optimoptions.s_weight = 10;
optimoptions.add_ccost = true; 

% ---- Initial Guess ----%
traj0 = load('TrajOpt_10-Oct-2017 21:09:05');

T = 100;
N = 21;
% x0 = hamr.getInitialState();
% qi = x0(1:nq); x1 = x0; x1(1) = 30; qf = x1(1:nq);
t_init = traj0.xtraj.getBreaks(); 
xi = traj0.xtraj.eval(t_init(1)); qi = xi(1:nq); 
xf = traj0.xtraj.eval(t_init(1)); qf = xf(1:nq); 

T_span = [0 2*T];
traj_opt = VariationalTrajectoryOptimization(hamr,N,T_span,optimoptions);

% load TrajOpt_10-Oct-2017 19:38:19

traj_init.x = traj0.xtraj;
traj_init.u = traj0.utraj;
traj_init.c = traj0.ctraj;
traj_init.b = traj0.btraj;
traj_init.psi = traj0.psitraj;
traj_init.eta =  traj0.etatraj;
% traj_init.s =  PPTrajectory(zoh(t_init(1:end-1),0.1*ones(1,N-1)));
% 
% traj_init.x = PPTrajectory(foh([0 T],[x0, x1]));
% traj_init.u = PPTrajectory(zoh(t_init(1:end-1),0.001*randn(nu,N-1)));
% traj_init.c = PPTrajectory(zoh(t_init(1:end-1),0.001*randn(traj_opt.nC,N-1)));
% traj_init.b = PPTrajectory(zoh(t_init(1:end-1),0.001*randn(traj_opt.nC*traj_opt.nD,N-1)));
% traj_init.psi = PPTrajectory(zoh(t_init(1:end-1),0.001*randn(traj_opt.nC,N-1)));
% traj_init.eta =  PPTrajectory(zoh(t_init(1:end-1),0.001*randn(traj_opt.nC*traj_opt.nD,N-1)));
% traj_init.s =  PPTrajectory(zoh(t_init(1:end-1),0.1*ones(1,N-1)));


% -- Costs ---%
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
% traj_opt = traj_opt.addFinalCost(@final_cost_fun);
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
    @periodic_constraint, cnstr_opts);

periodic_cnst_vel = periodic_cnst_vel.setName('periodicity_vel');

periodic_inds = {[traj_opt.h_inds(1); traj_opt.x_inds(:,1); traj_opt.x_inds(:,2); ...
    traj_opt.u_inds(:,1); traj_opt.c_inds(:,1); traj_opt.b_inds(:,1); traj_opt.jl_inds(:,1); ...
    traj_opt.kl_inds(:,1); traj_opt.h_inds(N-1); traj_opt.x_inds(:,N-1); traj_opt.x_inds(:,N); ...
    traj_opt.u_inds(:,N-1); traj_opt.kl_inds(:,N-1)]};

% ----- Bounding box on Height --- 
zlb = qi(3) - 1.5;
zub = qi(3) + 1.5;

% Add Initial/Final State Constraints 
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(qi(1:2)),1, 1:2);
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(zlb,zub),1:N, 3);
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(qf(1), Inf),N,1);

% Add Periodic Constraints
traj_opt = traj_opt.addPositionConstraint(periodic_cnst_pos,{[1 N]}, 2:nq);
traj_opt = traj_opt.addConstraint(periodic_cnst_vel, periodic_inds);

% Input constraints
traj_opt = traj_opt.addInputConstraint(BoundingBoxConstraint(umin, umax),1:N-1);


% Solver options
traj_opt = traj_opt.setSolver('snopt');
traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',10000);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',5000000);
traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',1000);
traj_opt = traj_opt.setSolverOptions('snopt','print','outputlog.txt');

disp('Solving...')
tic
[xtraj,utraj,ctraj,btraj,psitraj,etatraj,jltraj,kltraj,straj ...
    ,z,F,info,infeasible_constraint_name] = solveTraj(traj_opt,t_init,traj_init);
toc
%
    function [f,df] = running_cost_fun(h,x,u)
        R = (1/ulim)^2*eye(nu);
        g = (1/2)*u'*R*u;
        f = h*g;
        df = [g, zeros(1,nx), h*u'*R];
    end

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
        hNm1 = xin(1+2*nQ+nU+nC+nD*nC+nJL+1);
        qNm1 = xin(1+2*nQ+nU+nC+nD*nC+nJL+1+(1:nQ));
        qN = xin(1+2*nQ+nU+nC+nD*nC+nJL+1+nQ+(1:nQ));
        uNm1 = xin(1+2*nQ+nU+nC+nD*nC+nJL+1+2*nQ+(1:nU));
        klNm1 = xin(1+2*nQ+nU+nC+nD*nC+nJL+1+2*nQ+nU+(1:nKL));
        
        [p0, dp0] = left_legendre_transform_fun(traj_opt,h0,q0,q1,u0,c0,b0,jl0,kl0);
        [pN, dpN] = right_legendre_transform_fun(traj_opt,hNm1,qNm1,qN,uNm1,klNm1);
        
        f = pN - p0;
        
        df = zeros(nQ, numel(xin));
        df(:, 1:1+2*nQ+nU+nC+nD*nC+nJL+nKL) = -dp0;
        df(:, 1+2*nQ+nU+nC+nD*nC+nJL+nKL+(1:1+2*nQ+nU+nKL)) = dpN;
        
        fprintf('Periodic const: %f \r', f); 
    end

end
