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
x0 = zeros(100,1); x0(3) = 12.04;
q0 = x0(1:nq);

T0 = 1e3;       % one seconds
N = 11;

% ---- Initial Guess ----%
q1 = [50 ;q0(2:end)];  % move one body length
x1 = [q1;zeros(nv,1)];
t_init = linspace(0,T0,N);

% traj_init.x = PPTrajectory(foh([0 T0/2 T0],[x0, xm, x1]));
traj_init.x = PPTrajectory(foh([0 T0],[x0, x1]));
traj_init.u = PPTrajectory(zoh(t_init,0.1*randn(nu,N)));
T_span = [0 T0];

options.sweight = 1e4; 
traj_opt = VariationalTrajectoryOptimization(hamr,N,T_span,options);
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q0),1);
% traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(qm),7);
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q1(1)),N, 1);

z_ub = q0(3) + 2;
z_lb = q0(3) - 3;
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(z_lb,z_ub),2:N-1, 3);

% traj_opt = traj_opt.addVelocityConstraint(ConstantConstraint(zeros(nv,1)),1);

Flim = 200; % mN
traj_opt = traj_opt.addInputConstraint(BoundingBoxConstraint(-Flim*ones(nu,1),Flim*ones(nu,1)),1:N-1);

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
