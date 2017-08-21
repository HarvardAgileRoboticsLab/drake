function [hamr,xtraj,utraj,ctraj,btraj,...
    psitraj,etatraj,jltraj, kltraj, straj] = runHamrVariational()

urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMRVariational_scaled.urdf');

% options
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;

% floating (gnd contact) or in air (not floating)
% ISFLOAT = false;
ISFLOAT = true;
if ISFLOAT
    options.floating = ISFLOAT;
    options.collision = ISFLOAT;
    x0 = zeros(56,1); x0(3) = 20;
    %     options.terrain = RigidBodyFlatTerrain();
    
else
    options.floating = ISFLOAT;
    options.collision = ISFLOAT;
    x0 = zeros(44, 1);
    options.terrain = [];
end

% Build robot + visualizer
hamr = HamrVariational(urdf,options);
% 
% 
% simtraj = hamr.simulate([0, tf], x0);
% simtraj_scaled = PPTrajectory(foh(simtraj.getBreaks()*1e-3, simtraj.eval(simtraj.getBreaks())));
% simtraj_scaled = simtraj_scaled.setOutputFrame(simtraj.getOutputFrame());

v = hamr.constructVisualizer();
v.inspector(x0); 
tf = 10;

% v.playback(simtraj_scaled, struct('slider', true))

%%

% N=11;
% traj_init.x = PPTrajectory(foh([0 tf],[x0, x0]));
% t_init = linspace(0,tf,N);
% 
% 
% options.s_weight = 10;
% options.terrain = RigidBodyFlatTerrain();
% options.floating = true;
% options.ignore_self_collisions = true;
% options.use_bullet = false;
% nq = hamr.getNumPositions;
% 
% traj_opt = VariationalTrajectoryOptimization(hamr,N,tf,options);
% traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(x0(1:nq)),1);
% traj_opt = traj_opt.addVelocityConstraint(ConstantConstraint(x0(nq+(1:nq))),1);
% 
% % traj_opt = traj_opt.setSolver('ipopt');
% traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',1000000);
% traj_opt = traj_opt.setSolverOptions('snopt','MajorOptimalityTolerance',1e-5);
% traj_opt = traj_opt.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-5);
% traj_opt = traj_opt.setSolverOptions('snopt','MinorFeasibilityTolerance',1e-5);
% traj_opt = traj_opt.setSolverOptions('snopt','constraint_err_tol',1e-5);
% 
% traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);
% 
% tic
% [xtraj,utraj,ctraj,btraj,psitraj,etatraj,jltraj, kltraj, straj] = traj_opt.solveTraj(t_init,traj_init);
% toc
% 
% %%
% t = simtraj.getBreaks();
% x = simtraj.eval(simtraj.getBreaks());
% 
% figure(1); clf;
% plot(t, x(3,:)); hold on;
% 
%     function displayTraj(h,x,u)
%         
%         ts = [0;cumsum(h)];
%         for i=1:length(ts)
%             v.drawWrapper(0,x(:,i));
%             pause(h(1));
%         end
%         
%     end
% end