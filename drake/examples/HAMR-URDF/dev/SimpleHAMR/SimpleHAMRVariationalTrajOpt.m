function [hamr,x,u,c,b, psi,eta,jl, kl, s, z,F,info, ...
    infeasible_constraint_name] = SimpleHAMRVariationalTrajOpt(load_path, params)

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
nu = hamr.getNumInputs();
nc = hamr.getNumContactPairs();

% --- Set Input limits ---
ulim = 0.26;                 % set max force
umin = -ulim*ones(nu,1);
umax = ulim*ones(nu, 1);

% --- Initialize TrajOpt---
optimoptions.s_weight = params.s_weight;
optimoptions.add_ccost = true;
optimoptions.s_max = Inf; 

% parameters
T = 200;  % ms
r_des = params.r_des;         % mm
tht_des = params.tht_des;     % rad
xi = hamr.getInitialState();
xf = xi;
xf([1, 2, 6]) = [r_des*cos(tht_des); r_des*sin(tht_des); tht_des];
N = params.N;

if isempty(load_path)
    t_init = linspace(0, T, N);
    traj_init.x = PPTrajectory(foh([0 T],[xi, xf]));
    traj_init.u = PPTrajectory(zoh(t_init,0.0625*randn(nu,N)));
    traj_init.c = PPTrajectory(zoh(t_init,0.001*randn(nc,N)));
    traj_init.b = PPTrajectory(zoh(t_init,0.001*randn(4*nc,N)));
    traj_init.psi = PPTrajectory(zoh(t_init,0.001*randn(nc,N)));
    traj_init.eta = PPTrajectory(zoh(t_init,0.001*randn(4*nc,N)));
else
    traj_init = rmfield(load(load_path), {'jl', 'kl'});
    [t_init, traj_init] = re_interp(traj_init, N);
    traj_init.s = 0*traj_init.s + 1;
    
end

T_span = [0.5*T 1.5*T];
traj_opt = VariationalTrajectoryOptimization(hamr,N,T_span,optimoptions);

% -- Costs ---%
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
% traj_opt = traj_opt.addRunningCost(@foot_height_fun);
traj_opt = traj_opt.addFinalCost(@final_cost_fun);
traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);

% Add Initial State Constraints
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(xi(1:nq)),1);
traj_opt = traj_opt.addVelocityConstraint(ConstantConstraint(xi(nq+(1:nv))),1);

% final state constraints
lbf = xf(1:6) - [Inf; 5; 2; pi/6; pi/6; pi/6];
ubf = xf(1:6) + [Inf; 5; 2; pi/6; pi/6; pi/6];
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(lbf, ubf),N,1:6);

% z-bounding box
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(xi(3) - 2, xi(3) + 2),2:N,3);
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(xi(5) - pi/6, xi(5) + pi/6),2:N,5);


% Add joint limit constraints
jlmin = hamr.joint_limit_min; jlmax = hamr.joint_limit_max;
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(jlmin,jlmax),1:N);

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

    function [f,df] = running_cost_fun(h,xk,uk)
        xin = [h; xk; uk];
        [f, df] = running_cost(xin);
        fprintf('Running Cost: %d \r', f)
        
%         df_fd = zeros(size(df));
%         dxin = 1e-6*eye(length(xin));
%         for k = 1:length(xin)
%             df_fd(:,k) = (running_cost(xin+dxin(:,k)) ...
%                 - running_cost(xin-dxin(:,k)))/2e-6;
%         end
%         
%         disp('Running cost derivative error:');
%         disp(max(abs(df_fd(:)-df(:))));
%         
    end


    function [f, df] = running_cost(xin)
        
        h = xin(1);
        xk = xin(1+(1:nq+nv));
        uk = xin(1+nq+nv+(1:nu));
        
        Qdiag = zeros(nq + nv, 1);
        Qdiag(1:6) = [0.0; 1.0; 1.0; 10.0; 10.0; 10.0];
        Q = diag(Qdiag);
        
        R = 0.1*(1/ulim)^2*eye(nu);
        g = (1/2)*(uk'*R*uk + (xk - xf)'*Q*(xk - xf));
        f = h*g;
        df = [g, h*(xk - xf)'*Q, h*uk'*R];
        
    end


    function [f,df] = final_cost_fun(tN,xN)
        
        xin = [tN; xN];
        [f, df] = final_cost(xin);
        fprintf('Final Cost: %d \r', f) 
        
%         df_fd = zeros(size(df));
%         dxin = 1e-6*eye(length(xin));
%         for k = 1:length(xin)
%             df_fd(:,k) = (final_cost(xin+dxin(:,k)) ...
%                 - final_cost(xin-dxin(:,k)))/2e-6;
%         end
%         
%         disp('Final cost derivative error:');
%         disp(max(abs(df_fd(:)-df(:))));
        
    end

    function [f, df] = final_cost(xin)        
        tN = xin(1);
        xN = xin(1+(1:nq+nv));
        
        Qfdiag = zeros(nq + nv, 1);
        Qfdiag(1:6) = [10.0; 10.0; 0.1; 10.0; 10.0; 10.0];                 
        Qfdiag([7,8,12]) = 10; 
        Qf = diag(Qfdiag);
        f = (1/2)*(xN-xf)'*Qf*(xN-xf);
        df = [0, (xN-xf)'*Qf];        
    end

%     function [f,df] = foot_height_fun(h,x,u)
%         q = x(1:nq);
%         [phi,~,~,~,~,~,~,~,n] = hamr.contactConstraints(q,false,struct('terrain_only',true));
%         phi0 = [.5;.5;.5;.5];
%         K = 10;
%         I = find(phi < phi0);
%         f = K*(phi(I) - phi0(I))'*(phi(I) - phi0(I));
%         % phi: 2x1
%         % n: 2xnq
%         df = [0 2*K*(phi(I)-phi0(I))'*n(I,:) zeros(1,nv+nu)];
%     end

    function displayTraj(h,x,u)
        disp('Displaying Trajectory...')
        h = h/1e3;
        ts = [0;cumsum(h)];
        for i=1:length(ts)
            v.drawWrapper(0,x(:,i));
            pause(5*h(1));
        end
        
    end

    function [t_init, traj_init] = re_interp(traj_init, N)
        tt = traj_init.x.getBreaks();
        N0 = numel(tt);
        if N0 == N
            t_init = tt;
        else
            t_init = linspace(tt(1), tt(end), N);
            names = fieldnames(traj_init);
            for i = 1:numel(names)
                disp(i)
                if ~strcmpi(names{i}, 'params')
                    vals = interp1(tt', traj_init.(names{i}).eval(tt)',t_init')';
                    traj_init.(names{i}) = PPTrajectory(foh(t_init, vals));
                end
            end
        end
        
    end



end
