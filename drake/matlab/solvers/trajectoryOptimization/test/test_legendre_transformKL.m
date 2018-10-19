function [tt, xx] = test_legendre_transformKL()
options.floating = false;
options.use_bullet = false;
file = fullfile(getDrakePath,'examples','SimpleFourBar', 'FourBar_JointLimits3DOF.urdf');
plant = RigidBodyManipulator(file, options);  %RigidBodyManipulator(file,options);
% v = plant.constructVisualizer(); 

nq = plant.getNumPositions();
nv = plant.getNumVelocities();
nx = nq+nv;
nu = plant.getNumInputs();

% Nominal Data
qi = [4*pi/5; pi/5; 4*pi/5];  
vi = 0*qi; 
xi = [qi; vi]; 
ui = 0;
% qi = xi(1:nq);
% vi = 0*qi;%xi(nq+(1:nv));

% v.inspector([xi; vi]); 
% 
% w = warning('off','Drake:DrakeSystem:ConstraintsNotEnforced');
% xtraj = plant.simulate([0 10]);
% warning(w);

% v.playback(xtraj); 

T0 = 2;
N = 21;

% traj_init.x = PPTrajectory(foh([0 T0/2 T0],[x0, xm, x1]));
t_init = linspace(0, T0, N);
traj_init.x = PPTrajectory(foh([0, T0], [xi, xi]));
traj_init.u = PPTrajectory(zoh([0, T0], [ui, ui]));

T_span = [T0-1 2*T0];

optimoptions.s_weight = 10;
optimoptions.joint_limit_collisions = true;  
optimoptions.integration_method = VariationalTrajectoryOptimization.EULER; 

% optimoptions.periodic = 1;
traj_opt = VariationalTrajectoryOptimization(plant,N,T_span,optimoptions);


% Position 
Qperiodic = [-eye(nq), eye(nq)]; 
periodic_cnst_pos = LinearConstraint(zeros(nq,1),zeros(nq,1),Qperiodic);
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


% Add Periodic Constraints
traj_opt = traj_opt.addPositionConstraint(periodic_cnst_pos,{[1 N]});
traj_opt = traj_opt.addConstraint(periodic_cnst_vel, periodic_inds);

traj_opt = traj_opt.setSolver('snopt');
traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',10000);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',1000000);
traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',1000);

% traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);

disp('Solving...')
tic
[xtraj,utraj,ctraj,btraj, psitraj,etatraj,jltraj, kltraj,straj, ...
    z,F,info,infeasible_constraint_name] = traj_opt.solveTraj(t_init,traj_init);
toc
%%

tt = xtraj.getBreaks();
h = tt(end)/N;
xx = xtraj.eval(tt);
qq = xx(1,:);
xx2 = xtraj.eval(tt + h/2);
vv = xx2(2,:);

uu = utraj.eval(tt + h/2);

figure(1); clf;
subplot(2,1,1); hold on;
plot(tt, rad2deg(qq));
plot(tt+h/2, rad2deg(vv));
subplot(2,1,2); hold on;
plot(tt+h/2, uu(1,:));

v = plant.constructVisualizer();
v.playback(xtraj, struct('slider', true));

    function [f,df] = periodic_constraint_fun(h0,q0,q1,u0,c0,b0,jl0,kl0, ...
            hNm1,qNm1,qN,uNm1,klNm1)
        
        xin = [h0;q0;q1;u0;c0;b0;jl0;kl0; ...
            hNm1;qNm1;qN;uNm1;klNm1];      
        [f,df] = periodic_constraint(xin);
        %         f
        
%         df_fd = zeros(size(df));
%         step = sqrt(eps(max(xin)));
%         dxin = step*eye(length(xin));
%         for k = 1:length(xin)
%             xin + dxin(:,k); 
%             df_fd(:,k) = (periodic_constraint(xin+dxin(:,k)) - ...
%                 periodic_constraint(xin-dxin(:,k)))/(2*step);
%         end
%         fprintf('Periodic const: %f \r', max(abs(f))); 
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
end