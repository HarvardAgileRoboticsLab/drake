function xtraj=contactImplicitBall(plant,N,x0,integration)
% tests that the contact implicit trajectory optimization can reproduce a
% simulation of the falling brick
if nargin<1
options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.use_bullet = false;

%plant = PlanarRigidBodyManipulator(fullfile(getDrakePath,'matlab','systems','plants','test','FallingBrickContactPoints.urdf'),options);
%x0 = [0; 1; .1; 0; 0; 0];
file = fullfile(getDrakePath,'matlab','systems','plants','test','ball.urdf');
plant = RigidBodyManipulator(file,options);
end
if nargin < 2
  N=10;
end
if nargin<3
q0 = [0
      0
      1.0000
      0.1
     -0.1
      0.0];
    
x0 = [q0;0*q0];
end

tf=2.0;
dt = tf/(N-1);

p = TimeSteppingRigidBodyManipulator(plant,dt);

tic
xtraj = p.simulate([0 tf],x0);
toc

% v = p.constructVisualizer();
% v.playback(xtraj,struct('slider',true));

% %[0;0;.8;0.05*randn(3,1);zeros(6,1)];
% visualize=true;
% 
% % plant_ts = TimeSteppingRigidBodyManipulator(plant,tf/(N-1));
% % w = warning('off','Drake:TimeSteppingRigidBodyManipulator:ResolvingLCP');
% % xtraj_ts = simulate(plant_ts,[0 tf],x0);
% % x0 = xtraj_ts.eval(0);
% % warning(w);
% if visualize
%   v = constructVisualizer(plant);
% end
% 
% options = struct();
% if nargin > 3
%     if strcmp(integration, 'midpoint')
%         options.integration_method = ContactImplicitTrajectoryOptimization.MIDPOINT;
%     elseif strcmp(integration, 'backward')
%         options.integration_method = ContactImplicitTrajectoryOptimization.BACKWARD_EULER;
%     elseif strcmp(integration, 'forward')
%         options.integration_method = ContactImplicitTrajectoryOptimization.FORWARD_EULER;
%     else
%         options.integration_method = ContactImplicitTrajectoryOptimization.MIXED;
%     end
% else
%     options.integration_method = ContactImplicitTrajectoryOptimization.MIXED;
% end
% 
% scale_sequence = [1;.01;0.0001];
% 
% for i=1:length(scale_sequence)
%   scale = scale_sequence(i);
% 
%   options.compl_slack = scale*.01;
%   options.lincompl_slack = scale*.001;
%   options.jlcompl_slack = scale*.01;
%   
%   prog = ContactImplicitTrajectoryOptimization(plant,N,tf,options);
%   prog = prog.setSolverOptions('snopt','MajorIterationsLimit',200);
%   prog = prog.setSolverOptions('snopt','MinorIterationsLimit',200000);
%   prog = prog.setSolverOptions('snopt','IterationsLimit',200000);
%   % prog = prog.setCheckGrad(true);
%   
% %   snprint('snopt.out');
%   
%   % initial conditions constraint
%   prog = addStateConstraint(prog,ConstantConstraint(x0),1);
%   
%   if i == 1,
%     traj_init.x = PPTrajectory(foh([0,tf],[x0,x0]));
%   else
%     traj_init.x = xtraj;
%     traj_init.l = ltraj;
%   end
%   tic
%   [xtraj,utraj,ltraj,~,z,F,info] = solveTraj(prog,tf,traj_init);
%   toc
% end
% 
% if visualize
%   v.playback(xtraj);
% end

% ts_plant = TimeSteppingRigidBodyManipulator(plant,0.0005,options);

% sim_traj = ts_plant.simulate([0,tf],x0);

% 
% ts = sim_traj.getBreaks();
% %tt = xtraj.getBreaks();
% xx = xtraj.eval(ts);
% xs = sim_traj.eval(ts);
% 
% figure(1);
% for i=1:12
%   subplot(2,6,i);
%   plot(ts,xx(i,:),'b');
%   hold on;
%   plot(ts,xs(i,:),'r');
%   hold off;
% end


% % check if the two simulations did the same thing:
% ts = getBreaks(xtraj_ts);
% valuecheck(ts,getBreaks(xtraj));
% xtraj_data = xtraj.eval(ts); 
% xtraj_ts_data = xtraj_ts.eval(ts);
% nq = plant.getNumPositions();
% nv = plant.getNumVelocities();
% valuecheck(xtraj_data(1:nq,:),xtraj_ts_data(1:nq,:),position_tol); % is there a correct tolerance here?
% valuecheck(xtraj_data(nq+(1:nv),:),xtraj_ts_data(nq+(1:nv),:),velocity_tol); % is there a correct tolerance here?
% 
% % TIMEOUT 750
