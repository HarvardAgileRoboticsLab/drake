function [xtrajf,utraj,z,F,info,infeasible_constraint_name] = SLTrajOpt()
close all; clc
%% Build Single Leg
sl_urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'dev', 'SimpleHAMR', ...
    'urdf','SL_scaled.urdf');

% options
% options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = false;

SL = SLRBM(sl_urdf, options);
% SLnc = SL.removeAllStateConstraints(); 

v = SL.constructVisualizer();

nq = SL.getNumPositions();
nv = SL.getNumVelocities();
nx = nq + nv;
nu = SL.getNumInputs();
nl = SL.nl;

%% Set up Forward Sim

T = 100;
N = 21;
t_init = linspace(0, T, N);
% q0 = [  -1.1469e-06
%     3.4464e-08
%     -2.1919e-06
%     -2.5376e-03
%     -1.4714e-03
%     1.5752e-04
%     -2.6951e-03
%     -2.1793e-06];

% x0 = [q0; 0*q0];
% fd = 0.002;
% ul = 0.1*sin(2*pi*fd*t_init);
% figure; clf; plot(t_init, ul);
% us = zeros(size(ul));
% utraj = zoh(t_init, [ul; us]);
% 
% SLTSRBM = TimeSteppingRigidBodyManipulator(SL, 0.05, options); 
x0 = zeros(nx, 1);
% utraj = PPTrajectory(zoh(t_init, zeros(nu, N)));
% utraj = utraj.setOutputFrame(SLTSRBM.getInputFrame);
% SLPassive = cascade(utraj, SLTSRBM);
% xtraj = SLPassive.simulate(t_init, x0);
% xx = xtraj.eval(xtraj.getBreaks());
% qq = xtraj.eval(xtraj.getBreaks()); qq = qq(1:SL.getNumPositions(), :);
% qtraj_scaled = PPTrajectory(foh(xtraj.getBreaks()*1e-3, qq));
% qtraj_scaled = qtraj_scaled.setOutputFrame(v.getInputFrame());
% v.playback(qtraj_scaled, struct('slider', true));

% figure(1); clf; hold on;
% act_dof = SL.getActuatedJoints();
% plot(xtraj.getBreaks(), xx(act_dof(1),:));
% plot(xtraj.getBreaks(), xx(act_dof(2),:));
% legend('Swing', 'Lift');


%% Set up Traj Opt
optimopt.integration_method = 3;                            % Forward Euler
optimopt.time_option = 1;                                   % all steps are constant
traj_opt = DirtranTrajectoryOptimization(SL, N, T, optimopt);
% traj_opt = DircolTrajectoryOptimization(SL, N, T, optimopt);


%input constraint
ulim = [zeros(nu-nl, 1); Inf(nl, 1)]; 
traj_opt = traj_opt.addInputConstraint(BoundingBoxConstraint(-ulim, ulim), 1:N);

% state constraints
[qmin, qmax] = SL.getJointLimits();
traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(qmin, qmax), 1:N, 1:nq);

TOL = 1e-4;
loop_const = SL.loop_const; 
for ind = 1:6 %numel(loop_const)
    loopi = loop_const{ind};
    lb = loopi.lb; ub = loopi.ub;
%     if ind == 1
%         TOL = [TOLi; TOLi; 10];
%     elseif ind == 2
%         TOL = [10; TOLi; TOLi];
%     elseif ind == 3
%         TOL = [TOLi; 10; TOLi];
%     end
        
%         pci = pci.setBounds
    loopi = loopi.setBounds(lb-TOL, ub+TOL);
%     name = cell(size(pci.name));
%     for k = 1:numel(pci.name)
%         name{k} = [pci.name{k}, '[',num2str(ind),']'];
%     end
%     pci = pci.setName(name);
    traj_opt = traj_opt.addStateConstraint(loopi, 2:N, 1:nq);
end

traj_opt = traj_opt.addStateConstraint(ConstantConstraint(x0), 1);
traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);

% % Solver options
traj_opt = traj_opt.setSolver('snopt');
traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',10000);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',6000000);
traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',1000);
% traj_opt = traj_opt.setSolverOptions('snopt','print','outputlog.txt');
% 
% traj_opt = traj_opt.setSolverOptions('snopt','MajorOptimalityTolerance',1e-3);
% traj_opt = traj_opt.setSolverOptions('snopt','MinorOptimalityTolerance',1e-3);
% traj_opt = traj_opt.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-3);
% traj_opt = traj_opt.setSolverOptions('snopt','MinorFeasibilityTolerance',1e-3);
% traj_opt = traj_opt.setSolverOptions('snopt','constraint_err_tol',1e-3);

%% Init and Solve

traj_init.x = PPTrajectory(zoh([0, T], [x0, x0]));

% traj_init.x = PPTrajectory(foh(xtraj.getBreaks(), xx));
traj_init.u = PPTrajectory(zoh(t_init, zeros(nu, N)));

% xfoot = xfoot.eval(t_init);

% xtraj = SL.simulate(t_init, x0);
% utraj = [];
tic; 
[xtrajf,~,z,F,info,infeasible_constraint_name] = traj_opt.solveTraj(T, traj_init);
toc

tf = xtrajf.getBreaks();
qq = xtrajf.eval(tf); qq = qq(1:SL.getNumPositions(), :);
qtraj_scaled = PPTrajectory(foh(tf*1e-3, qq));
qtraj_scaled = qtraj_scaled.setOutputFrame(v.getInputFrame());



%%
xf = xtrajf.eval(t_init);
for j = 1:nq
    figure(j); clf; hold on;
%     plot(xtraj.getBreaks(), xx(j,:))
    plot(t_init, xf(j,:))
    title(SL.body(j+1).jointname)
end
tilefigs;

v.playback(qtraj_scaled, struct('slider', true));

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

