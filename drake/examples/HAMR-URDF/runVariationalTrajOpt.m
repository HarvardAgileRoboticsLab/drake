function [hamr,xtraj,utraj,ctraj,btraj,...
    psitraj,etatraj,jltraj, kltraj, straj] = runVariationalTrajOpt()

% file
urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMRVariational_scaled.urdf');

% options
options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = true;
options.collision = true;
% options.dt = 1; 

% ts_hamr = HamrVariationalTSRBM(urdf, options); 
hamr = HamrVariationalRBM(urdf,options);
v = hamr.constructVisualizer();

nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nx = nq+nv;
nu = hamr.getNumInputs();

% Nominal Data
nom_traj = load('NomTraj.mat');
% x0 = zeros(100,1); x0(3) = 12.04;
% q0 = x0(1:nq);

T0 = nom_traj.tt(end);
N = 11;

% ---- Initial Guess ----%
% q1 = [50 ;q0(2:end)];  % move one body length
% x1 = [q1;zeros(nv,1)];
t_init = linspace(0,T0,N);
x_init = interp1(nom_traj.tt, nom_traj.yy(1:2*nq, :)', t_init)'; 
u_init = interp1(nom_traj.tt, nom_traj.yy(2*nq+1:2*nq+nu, :)', t_init)'; 

% traj_init.x = PPTrajectory(foh([0 T0/2 T0],[x0, xm, x1]));
traj_init.x = PPTrajectory(foh(t_init, x_init));
traj_init.u = PPTrajectory(zoh(t_init, u_init));

T_span = [0 T0];

options.sweight = 10; 
traj_opt = VariationalTrajectoryOptimization(hamr,N,T_span,options);



traj_init.c = PPTrajectory(zoh([0, T0], [1e-3*rand(traj_opt.nC, 1), 1e-3*rand(traj_opt.nC, 1)])); 
traj_init.b = PPTrajectory(zoh([0, T0], [1e-3*rand(traj_opt.nD*traj_opt.nC, 1), 1e-3*rand(traj_opt.nD*traj_opt.nC, 1)])); 
traj_init.psi = PPTrajectory(zoh([0, T0], [1e-3*rand(traj_opt.nC, 1), 1e-3*rand(traj_opt.nC, 1)])); 
traj_init.eta = PPTrajectory(zoh([0, T0], [1e-3*rand(traj_opt.nC*traj_opt.nD, 1), 1e-3*rand(traj_opt.nC*traj_opt.nD, 1)])); 
traj_init.jl = PPTrajectory(zoh([0, T0], [1e-3*rand(traj_opt.nJL, 1), 1e-3*rand(traj_opt.nJL, 1)])); 
traj_init.kl = PPTrajectory(zoh([0, T0], [1e-3*rand(traj_opt.nKL, 1), 1e-3*rand(traj_opt.nKL, 1)])); 





% traj_opt = traj_opt.addRunningCost(@running_cost_fun);
% traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(xinit(),1);
% traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(qm),7);
% traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q1(1)),N, 1);
% 
% y_ub = q0(2) + 10;
% y_lb = q0(2) - 10;
% z_ub = q0(3) + 2;
% z_lb = q0(3) - 3;

% traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(...
%     [y_lb; z_lb] ,[z_ub; y_ub]),2:N-1, 2:3);

% traj_opt = traj_opt.addVelocityConstraint(ConstantConstraint(zeros(nv,1)),1);

% Flim = 0.2; % mN
% traj_opt = traj_opt.addInputConstraint(BoundingBoxConstraint(-Flim*ones(nu,1),Flim*ones(nu,1)),1:N-1);

traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',10000);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',1000000);
traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',1000);

traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);

disp('Solving...')
tic
[xtraj,utraj,ctraj,btraj, psitraj,etatraj,jltraj, kltraj,straj, ...
    z,F,info,infeasible_constraint_name] = traj_opt.solveTraj(t_init,traj_init);
toc


    function [f,df] = running_cost_fun(h,x,u)
        R = 2*(1/Flim)^2*eye(nu);
        g = (1/2)*u'*R*u;
        f = h*g;
        df = [g, zeros(1,nx), h*u'*R];
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
