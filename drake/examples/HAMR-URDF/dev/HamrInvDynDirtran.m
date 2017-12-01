function [hamr, xtraj,utraj,z,F,info] = HamrInvDynDirtran(tt, urdf, feet, feet_pos, td, options)
%% Build Single Leg

urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'urdf',[urdf, '.urdf']);

options.pf = feet_pos;
options.feet = feet;

hamr = HamrRBMID(urdf, options);
v = hamr.constructVisualizer();

nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nx = nq + nv;
nu = hamr.getNumInputs();
nl = hamr.nl;
nc = hamr.nc;
nd = hamr.nd;

%% Set up Trajectory Optimization

% -- IC --- %
T = tt(end);
No = numel(tt);
UPSAMPLE = 1;
N = UPSAMPLE*(No-1) + 1;
t_init = linspace(0, T, N);

% --- Init Dirtran --- %
optimopt.integration_method = 3;                            % Midpoint Euler
optimopt.time_option = 1;                                   % all steps are constant
traj_opt = DirtranTrajectoryOptimization(hamr, N, T, optimopt);

% --- objective ---- %
% traj_opt = traj_opt.addRunningCost(@running_cost_fun);

% ---- contact forces ---- %
cc = td.ctraj.eval(t_init);
cc = cc(:, 1:N-1);

bb = td.btraj.eval(t_init);
bb = bb(:,1:N-1);

% --- COM Positon and Velocity ---%

if options.floating
    xx = td.xtraj.eval(t_init);
    qq = xx(1:6,:);
    
    vv = xx(14+(1:6),:);
    vv = vv(:, 1:N-1);
end

% x0 = [xx(1:6,1); zeros(       % initial condition

% --- Shift quantities applied at middle of time step in Variational
Aperm = spdiags(0.5*ones(N, 2), 0:1, N-1, N);
cshift = (Aperm\cc')';
bshift = (Aperm\bb')';

if options.floating
    vshift = (Aperm\vv')';
end

ulim = 0.3;

% Fix foot and COM trajectories
XFEET = td.xfoot_traj.eval(t_init);


for ind = 1:N
    if ind < N
        umini = [-ulim*ones(nu-nl-nc*(1+nd),1); -Inf(nl,1); cshift(:,ind); bshift(:,ind)];
        umaxi = [ulim*ones(nu-nl-nc*(1+nd),1); Inf(nl,1); cshift(:,ind); bshift(:,ind)];
    else
        umini = [zeros(nu-nl-nc*(1+nd),1); zeros(nl,1); cshift(:,ind); bshift(:,ind)];
        umaxi = [zeros(nu-nl-nc*(1+nd),1); zeros(nl,1); cshift(:,ind); bshift(:,ind)];
    end
    traj_opt = traj_opt.addInputConstraint(BoundingBoxConstraint(umini, umaxi), ind);
    
    TxyzCOM = 5;
    TrpyCOM = deg2rad(50); 
    TCOM = [TxyzCOM*ones(3,1); TrpyCOM*ones(3,1)];
    
    TLeg = 0.05; 

    if ind > 2
        if options.floating
            traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(qq(:,ind)-TCOM, qq(:,ind)+TCOM), ...
            ind, (1:6));
        end
    
        foot_constrainti = FunctionHandleConstraint(-TLeg*ones(nc*size(feet_pos,1),1), ...
            TLeg*ones(nc*size(feet_pos,1),1), nq, @(x)foot_const_fun(hamr, x, XFEET(:,ind)));
        foot_constrainti = foot_constrainti.setName(sprintf('foot_constraint[%d]', ind)); 
        traj_opt = traj_opt.addStateConstraint(foot_constrainti, ind, 1:nq); 
    end    
%     traj_opt = traj_opt.addCost(FunctionHandleObjective(nu,@(x)objective_fun(x, qq(:,ind)),1), ...
%         traj_opt.u_inds(:,ind));
end


% joint limit constraint
[qmin, qmax] = hamr.getJointLimits();
traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(qmin, qmax), 2:N, 1:nq);

% loop const
loop_constraint = FunctionHandleConstraint(zeros(nl,1), zeros(nl,1), nq, @(x)loop_const(hamr, x));
traj_opt = traj_opt.addStateConstraint(loop_constraint, 2:N, 1:nq);

% Display function
traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);

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

if options.floating
    x0 = [qq(:,1); zeros(nq-6, 1); zeros(nv,1)];    
    x1 = [qq(:,N); zeros(nq-6, 1); zeros(nv,1)];
    traj_opt = traj_opt.addStateConstraint(ConstantConstraint(x0), 1); 
    traj_init.x = PPTrajectory(foh([0, T], [x0, x1]));
else
    x0 = zeros(nx,1); 
    traj_opt = traj_opt.addStateConstraint(ConstantConstraint(x0), 1);
    traj_init.x = PPTrajectory(zoh([0, T], [x0, x0]));

end

traj_init.u = PPTrajectory(zoh(t_init, 0.01*rand(nu, N)));

[xtraj,utraj,z,F,info] = traj_opt.solveTraj(T, traj_init);

%% Costs and Constraints
% 
%     function [f,df] = running_cost_fun(h,x,u)
%         
%         R = (1/ulim)^2*eye(nu-nc*(1+nd));
%         g = (1/2)*u(1:(nu-nc*(1+nd)))'*R*u(1:(nu-nc*(1+nd)));
%         f = h*g;
%         df = [g, zeros(1,nx), h*u(1:nu-nc*(1+nd))'*R, zeros(1,nc*(1+nd))];
%     end

    function [f, df] = objective_fun(u, qq)
%         u0 = [zeros(nu-nl-nc*(1+nd), 1); zeros(nl, 1); cc; bb];
%         R = blkdiag(eye(nu-nl-nc*(1+nd)), zeros(nl), 10*eye(nc), 10*eye(nc*nd));   
%         f = (1/2)*(u-u0)'*R*(u-u0); 
%         df = (u-u0)'*R; 
    end

    function [f, df] = foot_const_fun(obj, xin, xfeetd)
        
        %
        [f,df] = foot_const(obj, xin, xfeetd);
        fprintf('Foot Constraint: %f \r', max(abs(f)));
        
%         df_fd = zeros(size(df));
%         step = sqrt(eps(max(xin)));
%         dxin = step*eye(length(xin));
%         for k = 1:length(xin)
%             df_fd(:,k) = (foot_const(obj, xin+dxin(:,k), xfeetd) - ...
%                 foot_const(obj, xin-dxin(:,k), xfeetd))/(2*step);
%         end
%         
% %         disp('Objective Derivative Error:');
%         disp(max(abs(df_fd(:)-df(:))));
        
    end

    function [f, df] =  foot_const(obj, x, xfeetd)
        
        PF = obj.pf; % position of foot in local frame        
        
            qk = x(1:nq);
            vk = zeros(nv, 1);
            
            kinsol = obj.doKinematics(qk, vk, struct('compute_gradients', true));
            
            xfooti = zeros(size(PF,1)*nc, 1);
            for j = 1:obj.nc
                [xfooti(size(PF,1)*(j-1)+1:size(PF,1)*j), Jk(size(PF,1)*(j-1)+1:size(PF,1)*j, :)] = ...
                    obj.forwardKin(kinsol, obj.findLinkId(obj.feet{j}), PF(:,j));
            end
            
            f = xfooti - xfeetd;      % position
            df = Jk; 
        
    end

    function [f, df] = loop_const(obj, q)
        [k, K] = obj.positionConstraints(q);
        f = k(obj.valid_loops);
        df = K(obj.valid_loops, :);
    end

    function displayTraj(h,x,u)
        h = h/1e3;
        ts = [0;cumsum(h)];
        for i=1:length(ts)
            v.drawWrapper(0,x(:,i));
            pause(5*h(1));
        end
        
    end
end

