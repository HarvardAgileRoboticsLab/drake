function [SL, xtraj,utraj,z,F,info] = SLInverseDynDirtran(tt, urdf, foot, pf, xfoot, ctraj, btraj)
%% Build Single Leg

sl_urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'dev', 'SimpleHAMR', ...
    'urdf',[urdf, '.urdf']);

% options
options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = false;
options.pf = pf;
options.foot = foot; 

SL = SLRBM(sl_urdf, options);
v = SL.constructVisualizer();

nq = SL.getNumPositions();
nv = SL.getNumVelocities();
nx = nq + nv;
nu = SL.getNumInputs();
nl = SL.nl;
nc = SL.nc;
nd = SL.nd;

%% Set up Trajectory Optimization
T = tt(end);
No = numel(tt);
UPSAMPLE = 2;
N = UPSAMPLE*(No-1) + 1;
t_init = linspace(0, T, N);
x0 = zeros(nx, 1);


optimopt.integration_method = 3;                            % Midpoint Euler
optimopt.time_option = 1;                                   % all steps are constant
traj_opt = DirtranTrajectoryOptimization(SL, N, T, optimopt);

% objective
traj_opt = traj_opt.addCost(FunctionHandleObjective(N*nx,@(x)objective_fun(SL,x),1), ...
    traj_opt.x_inds(:));
% traj_opt = traj_opt.addRunningCost(@running_cost_fun);

% contact forces
cc = ctraj.eval(t_init);
cc = cc(:, 1:N-1);


bb = btraj.eval(t_init);
% bb = interp1(tt', bb', t_init')';
% bb = reshape(repmat(bb, UPSAMPLE,1), nc*nd, []);
bb = bb(:,1:N-1);

Aperm = spdiags(0.5*ones(N, 2), 0:1, N-1, N);
uc = (Aperm\cc')'; 
% uc = interp1(tt', uc', t_init')';

ub = (Aperm\bb')';
% ub = interp1(tt', ub', t_init')';

ulim = 0.4;

for ind = 1:N
    if ind < N
        umini = [-ulim*ones(nu-nl-nc*(1+nd),1); -Inf(nl,1); uc(:,ind); ub(:,ind)];
        umaxi = [ulim*ones(nu-nl-nc*(1+nd),1); Inf(nl,1); uc(:,ind); ub(:,ind)];        
    else
        umini = [zeros(nu-nl-nc*(1+nd),1); -Inf(nl,1); uc(:,ind); ub(:,ind)];
        umaxi = [zeros(nu-nl-nc*(1+nd),1); Inf(nl,1); uc(:,ind); ub(:,ind)];    
    end
    traj_opt = traj_opt.addInputConstraint(BoundingBoxConstraint(umini, umaxi), ind);
end

% joint limit constraint
[qmin, qmax] = SL.getJointLimits();
traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(qmin, qmax), 2:N, 1:nq);
traj_opt = traj_opt.addStateConstraint(ConstantConstraint(x0), 1);

% loop const
TOL = 1e-4;
loop_const = SL.loop_const;
for ind = 1:numel(loop_const)
    loopi = loop_const{ind};
    lb = loopi.lb; ub = loopi.ub;
    loopi = loopi.setBounds(lb-TOL, ub+TOL);
    traj_opt = traj_opt.addStateConstraint(loopi, 2:N, 1:nq);
end

% Display function
% traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);

% Solver options
traj_opt = traj_opt.setSolver('snopt');
traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',10000);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',5000000);
traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',1000);

traj_opt = traj_opt.setSolverOptions('snopt','MajorOptimalityTolerance',1e-4);
traj_opt = traj_opt.setSolverOptions('snopt','MinorOptimalityTolerance',1e-4);
traj_opt = traj_opt.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-4);
traj_opt = traj_opt.setSolverOptions('snopt','MinorFeasibilityTolerance',1e-4);
traj_opt = traj_opt.setSolverOptions('snopt','constraint_err_tol',1e-4);

%% Init and Solve

traj_init.x = PPTrajectory(zoh([0, T], [x0, x0]));
traj_init.u = PPTrajectory(zoh(t_init, zeros(nu,N)));

xfoot = xfoot.eval(t_init);

[xtraj,utraj,z,F,info] = traj_opt.solveTraj(T, traj_init);

%% Costs and Constraints

    function [f,df] = running_cost_fun(h,x,u)
        
        R = (1/ulim)^2*eye(nu-nl-nc*(1+nd));
        g = (1/2)*u(1:(nu-nl-nc*(1+nd)))'*R*u(1:(nu-nl-nc*(1+nd)));
        f = h*g;
        df = [g, zeros(1,nx), h*u(1:nu-nl-nc*(1+nd))'*R, zeros(1,nl+nc*(1+nd))];
    end

    function [f, df] = objective_fun(obj, xin)
        
        %
        [f,df] = objective(obj,xin);
        fprintf('Objective: %f \r', f);
        %
        %         df_fd = zeros(size(df));
        %         step = sqrt(eps(max(xin)));
        %         dxin = step*eye(length(xin));
        %         for k = 1:length(xin)
        %             df_fd(:,k) = (objective(obj, xin+dxin(:,k)) - objective(obj, xin-dxin(:,k)))/(2*step);
        %         end
        %
        %         disp('Objective Derivative Error:');
        %         disp(max(abs(df_fd(:)-df(:))));
        %
    end

    function [f, df] = objective(obj, x)
        
        PF = obj.pf; % position of foot in local frame
        
        a = 1; %1/(norm(PF)^2)/N;
        f = 0;
        df = zeros(1, numel(x));
        
        for i = 1:N
            qk = x(nx*(i-1)+(1:nq));
            vk = zeros(nv, 1);
            
            kinsoli = obj.doKinematics(qk, vk, struct('compute_gradients', true));
            [xfooti, Jk] = obj.forwardKin(kinsoli, obj.findLinkId(obj.foot), PF);
            
            f = f + (a/2)*(xfooti - xfoot(:,i))'*(xfooti - xfoot(:,i));      % position
            df(nx*(i-1)+(1:nx)) = [a*(xfooti - xfoot(:,i))'*Jk, zeros(1,nv)];      % df/dq, df/dv
        end
        
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

