function [hamr,xtraj,utraj,ctraj,btraj,...
    psitraj,etatraj,jltraj, kltraj, straj, ...
    z,F,info,infeasible_constraint_name] = HAMRVariationalTrajOpt(fname, phid)

% file
urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2.urdf');

% options
options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = true;
options.collision = true;

hamr = HamrRBM(urdf,options);
v = hamr.constructVisualizer();

% state/input dimenisons
nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nx = nq+nv;
nu = hamr.getNumInputs();
nl = hamr.nl; 

% --- Set Input limits ---
FlimS = 0.30;                 % set max force
FLimL = 0.30;
umin = -reshape([FlimS; FLimL]*ones(1, nu/2), nu, 1);
umax =  reshape([FlimS; FLimL]*ones(1, nu/2), nu, 1);
hamr = hamr.setInputLimits(umin, umax);

% --- Initialize TrajOpt---
optimoptions.s_weight = 200;
optimoptions.joint_limit_collisions = false;

% ---- Initial Guess ----%
if ~ isempty(phiS)
trajCOM = load(fname);                  % center of mass + contact forces
trajTrans = load([fname, '_fullRobot']);    % transmission + input + joint limits
trajVar = load([fname, '_Variational']);    % full planned traj

% Time 
% tt = trajTrans.xtraj.FL_scaled.getBreaks();
% tt = trajTrans.xtraj.getBreaks(); 
% hh = mean(diff(tt)); 
tt = trajVar.xtraj.getBreaks(); 
N = numel(tt);
tt = tt(1:N); 
T = tt(N);
T_span = [T T];

% % Build x-trajectory 
% xxCOM = trajCOM.xtraj.eval(tt); 
% xxFL = trajTrans.xtraj.FL_scaled.eval(tt);
% xxRL = trajTrans.xtraj.RL_scaled.eval(tt);
% xxFR = trajTrans.xtraj.FR_scaled.eval(tt);
% xxRR = trajTrans.xtraj.RR_scaled.eval(tt);
% 
xx = trajVar.xtraj.eval(tt); 
% xx = [xxCOM(1:6,:); xxFL(1:8,:); xxRL(1:8,:); xxFR(1:8,:); xxRR(1:8,:);
%     xxCOM(14+(1:6),:); xxFL(8+(1:8),:); xxRL(8+(1:8),:); xxFR(8+(1:8),:); xxRR(8+(1:8),:)];    
% 
% % Build u-trajectory 
% utraj = trajTrans.utraj;
% nut = size(utraj.FL_scaled.eval(tt), 1);
% uuTrans = [utraj.FL_scaled.eval(tt+hh/2);
%     utraj.RL_scaled.eval(tt+hh/2);
%     utraj.FR_scaled.eval(tt+hh/2);
%     utraj.RR_scaled.eval(tt+hh/2)];
% uu = [uuTrans(1:2, :); uuTrans(nut+(1:2),:); uuTrans(2*nut+(1:2),:); uuTrans(3*nut+(1:2), :)];
% 
% % Build loop force trajectory 
% ll = [uuTrans(2+(1:nl/4), :); uuTrans(nut+2+(1:nl/4),:); ...
%     uuTrans(2*nut+2+(1:nl/4),:); uuTrans(3*nut+2+(1:nl/4), :)];

% Initialzie optimization 
traj_opt = VariationalTrajectoryOptimization(hamr,N,T_span,optimoptions);
% traj_init.x = trajTrans.xtraj; xx = trajTrans.xtraj.eval(tt); 
% traj_init.u = trajTrans.utraj;
% traj_init.c = trajTrans.ctraj;
% traj_init.b = trajTrans.btraj;
% traj_init.psi= trajTrans.psitraj; 
% traj_init.eta = trajTrans.etatraj;
% traj_init.kl = trajTrans.kltraj;
% traj_init.x = PPTrajectory(foh(tt, xx)); 
% traj_init.u = PPTrajectory(zoh(tt, uu));
% traj_init.c = trajCOM.ctraj;
% traj_init.b = trajCOM.btraj;
% traj_init.psi = trajCOM.psitraj;
% traj_init.eta =  trajCOM.etatraj;
% traj_init.kl =  PPTrajectory(zoh(tt, ll));
% traj_init.s = PPTrajectory(zoh([0, T], [1, 1]));
traj_init.x = trajVar.xtraj; 
traj_init.u = trajVar.utraj;
traj_init.c = trajVar.ctraj;
traj_init.b = trajVar.btraj;
traj_init.psi = trajVar.psitraj;
traj_init.eta =  trajVar.etatraj;
traj_init.kl =  trajVar.kltraj; 
traj_init.s = PPTrajectory(zoh([0, T], [1, 1]));
else
    
    
end


% -- Costs ---%
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
for ind = 1:N
    traj_opt = traj_opt.addCost(FunctionHandleObjective(nq,@(x)tracking_cost(x, xx(1:nq,ind), phid(:,ind)),1), ...
        traj_opt.x_inds(1:nq,ind));
end
% traj_opt = traj_opt.addFinalCost(@final_cost_fun);

% -- Constraints ---%
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(xx(1:nq, 1)),1);
traj_opt = traj_opt.addVelocityConstraint(ConstantConstraint(0*xx(1:nq, 1)),1);
traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);

jlmin = hamr.joint_limit_min; %jlmin(3) = q0(3) - 1.5;
jlmax = hamr.joint_limit_max; %jlmax(3) = q0(3) + 1.5;

traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(jlmin,jlmax),2:N);
traj_opt = traj_opt.addInputConstraint(BoundingBoxConstraint(umin, umax),1:N-1);

% Solver options
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
    z,F,info,infeasible_constraint_name] = traj_opt.solveTraj(tt,traj_init);
toc

    function [f,df] = running_cost_fun(h,x,u)
        R = (2/(FlimS+FLimL))^2*eye(nu);
        g = (1/2)*u'*R*u;
        f = h*g;
        df = [g, zeros(1,nx), h*u'*R];
    end

    function [f, df] = tracking_cost(x, xd, phid)
        dim = 3;            % RigidBodyPose
        
        Q1 = eye(dim);
        Q2 = eye(numel(phid)); 
        
        kinsol = doKinematics(hamr, x, 0*x); 
        [phi,~,~,~,~,~,~,~,n] = hamr.contactConstraints(kinsol);         
        
        f = (1/2)*(x(1:dim)-xd(1:dim))'*Q1*(x(1:dim)-xd(1:dim)) + (1/2)*(phi-phid)'*Q2*(phi-phid);
        df = [(x(1:dim)-xd(1:dim))'*Q1, zeros(1, nq-dim)] + (phi-phid)'*Q2*n;  
        fprintf('Objective: %f \r', f)
    end
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
