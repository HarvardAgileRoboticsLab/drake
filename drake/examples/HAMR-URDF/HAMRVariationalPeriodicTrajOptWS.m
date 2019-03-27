function [hamr,xtraj,utraj,ctraj,btraj,...
    psitraj,etatraj,jltraj, kltraj, straj, ...
    z,F,info,infeasible_constraint_name] = HAMRVariationalPeriodicTrajOptWS(traj_init)

% file
urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2_TYM.urdf');

% options
options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = true; %false;
options.collision = true; %false;

hamr = HamrRBM(urdf,options);
v = hamr.constructVisualizer();

% state/input dimenisons
nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nx = nq+nv;
nu = hamr.getNumInputs();

% --- Set Input limits ---
FlimS = 0.28;                 % set max force
FLimL = 0.36;
umin = -reshape([FlimS; FLimL]*ones(1, nu/2), nu, 1);
umax =  reshape([FlimS; FLimL]*ones(1, nu/2), nu, 1);
hamr = hamr.setInputLimits(umin, umax);

% --- Initialize TrajOpt---
optimoptions.s_weight = 50;
optimoptions.joint_limit_collisions = false;
optimoptions.add_ccost = true;

% ---- Initial Guess ----%
tt = traj_init.t;
T = tt(end);
N = 51;
t_init = linspace(0, T, N);
xx0 = traj_init.x.eval(t_init(1));
xxf = traj_init.x.eval(t_init(N)); 


t_init = linspace(0, T, N);

% traj_init.x = traj_init.x.eval(t_init);
% traj_init.u = traj_init.u.eval(t_init);
% traj_init.c = traj_init.c.eval(t_init);
% traj_init.b = traj_init.b.eval(t_init);
% traj_init.psi = traj_init.psi.eval(t_init);
% traj_init.eta = traj_init.eta.eval(t_init);
% traj_init.kl = traj_init.kl.eval(t_init);

T_span = [0.9*T 1.1*T];
traj_opt = VariationalTrajectoryOptimization(hamr,N,T_span,optimoptions);

% -- Costs ---%
traj_opt = traj_opt.ad4dFinalCost(@final_cost);                          % close to final goal
% traj_opt = traj_opt.addRunningCost(@running_cost); 
traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);         % display

% 
% xxd = traj_init.x.eval(t_init); 
% % track nominal trajectory (loose)
% for ind = 1:N
%     traj_opt = traj_opt.addCost(FunctionHandleObjective(nq,@(x)tracking_cost(x, xxd(1:nq,ind)),1), ...
%         traj_opt.x_inds(1:nq,ind));
% end

% % try to get smooth u
% for ind = 1:N-2
%     traj_opt = traj_opt.addCost(FunctionHandleObjective(2*nu,@(u1, u2)udiff_cost(u1, u2)), ...
%         {traj_opt.u_inds(:,ind); traj_opt.u_inds(:,ind+1)});
% end

% % try to get smooth accelerations
% for ind = 2:N-1
%     traj_opt = traj_opt.addCost(FunctionHandleObjective(2*nv,@vdiff_cost_fun), ...
%         {traj_opt.h_inds(ind-1); traj_opt.h_inds(ind); traj_opt.x_inds(:,ind-1); ...
%         traj_opt.x_inds(:,ind); traj_opt.x_inds(:,ind+1)});
% end


%  `
% nqact =  numel(hamr.getActuatedJoints);
% nvarsCoT = nqact *(N-1) + nu*(N-1) +1;
% q_actuated = traj_opt.x_inds(hamr.getActuatedJoints,2:N);
% CoT_inds = {q_actuated(:); traj_opt.u_inds(:); traj_opt.x_inds(1,N)};
%
% traj_opt = traj_opt.addCost(FunctionHandleObjective(nvarsCoT,...
%     @(qact, u, xf)CoT_cost_fun(qact, u, xf)),CoT_inds);

% -----Periodic Constraint-----

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

% Add Initial/Final State Constraints
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(xx0(1:2)),1, 1:2);
% traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(xxf(1)-1, Inf),N,1);

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


    function [c,dc] = final_cost(tf,x)
        Q = eye(nq+nv); Q(1,1) = 1; 
        c = (1/2)*(x-xxf)'*Q'*(x-xxf);
        dc = [0, (x-xxf)'*Q];
        fprintf('Final Cost %f \r', c); 
    end

    function [c, dc] = vdiff_cost_fun(h1, h2, q1, q2, q3)
        
        xin = [h1; h2; q1; q2; q3];
        [c,dc] = vdiff_cost(xin);
%         
%         dc_fd = zeros(size(dc));
%         step = 1e-7;
%         dxin = step*eye(length(xin));
%         for k = 1:length(xin)
%             dc_fd(:,k) = (vdiff_cost(xin+dxin(:,k)) - vdiff_cost(xin-dxin(:,k)))/(2*step);
%         end
%         disp('VDiff Derivative Error:');
%         disp(max(abs(dc_fd(:)-dc(:))));
%         
        
    end

    function [c, dc] = vdiff_cost(xin)
        
        h1 = xin(1);
        h2 = xin(2);
        q1 = xin(2+(1:nq));
        q2 = xin(2+nq+(1:nq));
        q3 = xin(2+2*nq+(1:nq));
        
        %Take care of angle wrap-around
        vm1 = traj_opt.qdiff(q1,q2,h1);
        vm2 = traj_opt.qdiff(q2,q3,h2);
        
        a = 5; 
        vdiff = (vm1 - vm2);
        c = 0.5*a*(vdiff'*vdiff);
        I = eye(nv);
        dc = a*vdiff'*[-vm1/h1, vm2/h2, -I/h1, I/h1+I/h2, -I/h2];
%         fprintf('Vdiff Cost: %f \r', c)
    end

%     function [c,dc] = running_cost(h, x, u)        
%         
%         R = 0.1*(1/FLimL)^2*eye(nu);                         
%         g =(1/2)*(u'*R*u);
%         c = h*g;
%         dc = [g, zeros(1, numel(x)), h*u'*R];
%     end


    function [c, dc] = tracking_cost(x, xd)
        Q = 0.1*eye(nq);
        c = (1/2)*(x-xd)'*Q*(x-xd);
        dc = (x-xd)'*Q;
%         fprintf('Tracking Cost: %f \r', c)
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

% 
    function [f, df] = udiff_cost(u1, u2)
        udiff = u1-u2;
        f = 0.5*(udiff'*udiff);
        I = eye(length(u1));
        df = [udiff'*I,-udiff'*I];
    end



%     function [f, df] = running_cost_fun(h, x, u)
%
%         xin = [h; x; u];
%         [f,df] = running_cost(xin);
%
% %         df_fd = zeros(size(df));
% %         step = sqrt(eps(max(xin)));
% %         dxin = step*eye(length(xin));
% %         for k = 1:length(xin)
% %             df_fd(:,k) = (running_cost(xin+dxin(:,k)) - running_cost(xin-dxin(:,k)))/(2*step);
% %         end
% %
% %         disp('Running cost derivative error:');
% %         disp(max(abs(df_fd(:)-df(:))));
%     end
%
%
%     function [f, df] = running_cost(xin)
%
%         h = xin(1);
%         q = xin(1+(1:nq));
%         qd = xin(1+nq+(1:nv));
%
%         % lift hind legs up
%         kinsol = hamr.doKinematics(q, qd);
%         [phi,~,~,~,~,~,~,~,n] = hamr.contactConstraints(kinsol);
%         a = -0.1;
%         lc = [0; a; 0; a];
%
%         % minimize input torque
%         R = eye(nu);
%
%         f = lc'*phi;
%         fprintf('Hind Leg Cost: %f \r', f);
%         df = [0, lc'*n, zeros(1,nv), zeros(1,nu)];
%
% %         f = lc'*phi + (1/2)*u'*R*u;
% %         df = [0, lc'*n, zeros(1,nv), u'*R];
%     end

end
