function [hamr,xtraj,utraj,ctraj,btraj,...
    psitraj,etatraj,jltraj, kltraj, straj,z,F,info,infeasible_constraint_name] = runVariationalTrajOpt(LOAD)

% file
urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMRVariational_scaledV2.urdf');

% options
options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = true;
options.collision = true;

hamr = HamrVariationalRBM(urdf,options);
hamr = hamr.setJointLimits(-Inf(hamr.getNumPositions(), 1), Inf(hamr.getNumPositions(), 1));
hamr = compile(hamr);
v = hamr.constructVisualizer();

% state/input dimenisons
nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nx = nq+nv;
nu = hamr.getNumInputs();

% --- Set Input limits ---
FlimS = 0.375;                 % set max force
FLimL = 0.375; 
umin = -reshape([FlimS; FLimL]*ones(1, nu/2), nu, 1); 
umax =  reshape([FlimS; FLimL]*ones(1, nu/2), nu, 1); 
hamr = hamr.setInputLimits(umin, umax); 

% --- Initialize TrajOpt---
optimoptions.s_weight = 0; %1e-5; 
% optimoptions.joint_limit_collisions = false; 

% ---- Initial Guess ----%
if LOAD == 0   
        
    T = 100;                      
    N = 11; 
    x0 = hamr.getInitialState(); 
    load q0; 
    x0(1:nq) = q0;  
    q0 = x0(1:nq); v0 = x0(nq+1:(nq+nv));
    x1 = x0; %x1(1) = 10;    
    t_init = linspace(0,T,N);  
    T_span = [T T];    
    traj_opt = VariationalTrajectoryOptimization(hamr,N,T_span,optimoptions);    
    
    traj_init.x = PPTrajectory(foh([0 T],[x0, x1]));
    traj_init.u = PPTrajectory(zoh(t_init,0.001*randn(nu,N)));
    traj_init.c = PPTrajectory(zoh(t_init,0.001*randn(traj_opt.nC,N))); 
    traj_init.b = PPTrajectory(zoh(t_init,0.001*randn(traj_opt.nC*traj_opt.nD,N))); 
    traj_init.psi = PPTrajectory(zoh(t_init,0.001*randn(traj_opt.nC,N))); 
    traj_init.eta =  PPTrajectory(zoh(t_init,0.001*randn(traj_opt.nC*traj_opt.nD,N))); 
%     traj_init.jl = PPTrajectory(zoh([0, T], [1e-3*rand(traj_opt.nJL, 1), 1e-3*rand(traj_opt.nJL, 1)])); 
    traj_init.kl =  PPTrajectory(zoh(t_init,0.001*randn(traj_opt.nKL,N))); 
    traj_init.s = PPTrajectory(zoh([0, T], [1, 1]));  
else    
    prev_traj = load('none_var_150V_1Hz_60'); 
    
    N = 11; 
    T = 100; %prev_traj.tt(end);                      
    x0 = prev_traj.xx(:,1); 
    q0 = x0(1:nq); v0 = x0(nq+1:end); 
    t_init = linspace(0,T,N);  
    T_span = [T T];    
    traj_opt = VariationalTrajectoryOptimization(hamr,N,T_span,optimoptions);    
    
    traj_init.x = PPTrajectory(foh(prev_traj.tt, prev_traj.xx));
    traj_init.u = PPTrajectory(zoh(prev_traj.tt, prev_traj.uu));
    traj_init.c = PPTrajectory(zoh(prev_traj.tt, prev_traj.c_traj)); 
    traj_init.b = PPTrajectory(zoh(prev_traj.tt, prev_traj.beta_traj)); 
    traj_init.psi = PPTrajectory(zoh(prev_traj.tt, prev_traj.psi_traj)); 
    traj_init.eta =  PPTrajectory(zoh(prev_traj.tt, prev_traj.eta_traj)); 
%     traj_init.jl = PPTrajectory(zoh([0, T], [1e-3*rand(traj_opt.nJL, 1), 1e-3*rand(traj_opt.nJL, 1)])); 
    traj_init.kl =  PPTrajectory(zoh(prev_traj.tt, prev_traj.kl_traj));     
    traj_init.s = PPTrajectory(zoh([0, T], [1, 1]));  
end


    

% -- Costs ---%
% traj_opt = traj_opt.addRunningCost(@running_cost_fun);
% traj_opt = traj_opt.addFinalCost(@final_cost_fun);

% -- Constraints ---%
% traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q0),1);
% traj_opt = traj_opt.addVelocityConstraint(ConstantConstraint(v0),1);
traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);

% jlmin = hamr.joint_limit_min; jlmin(3) = q0(3) - 1.5;
% jlmax = hamr.joint_limit_max; jlmax(3) = q0(3) + 1.5;

% traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(jlmin,jlmax),2:N);
% traj_opt = traj_opt.addInputConstraint(BoundingBoxConstraint(umin, umax),1:N-1);

% Solver options
traj_opt = traj_opt.setSolver('snopt');
traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',10000);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',1000000);
traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',1000);
traj_opt = traj_opt.setSolverOptions('snopt','print','outputlog.txt');

traj_opt = traj_opt.setSolverOptions('snopt','MajorOptimalityTolerance',1e-3);
traj_opt = traj_opt.setSolverOptions('snopt','MinorOptimalityTolerance',1e-3);
traj_opt = traj_opt.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-3);
traj_opt = traj_opt.setSolverOptions('snopt','MinorFeasibilityTolerance',1e-3);
traj_opt = traj_opt.setSolverOptions('snopt','constraint_err_tol',1e-3);


disp('Solving...')
tic
[xtraj,utraj,ctraj,btraj, psitraj,etatraj,jltraj, kltraj,straj, ...
    z,F,info,infeasible_constraint_name] = traj_opt.solveTraj(t_init,traj_init);
toc

%     function [f,df] = running_cost_fun(h,x,u)
%         R = (2/(FlimS+FLimL))^2*eye(nu);
%         g = (1/2)*u'*R*u;
%         f = h*g;
%         df = [g, zeros(1,nx), h*u'*R];
%     end
% 
% 
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
end
