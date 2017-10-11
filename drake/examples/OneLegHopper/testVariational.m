warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
warning('off','Drake:RigidBodyManipulator:WeldedLinkInd');
warning('off','Drake:RigidBodyManipulator:UnsupportedJointLimits');
options.twoD = true;
%options.view = 'right';
options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.ignore_self_collisions = true;
options.use_bullet = false;
s = 'OneLegHopper.urdf';
m = PlanarRigidBodyManipulator(s,options);
dt = 0.0005;
r = TimeSteppingRigidBodyManipulator(m,dt,options);
v = m.constructVisualizer();

nQ = r.getNumPositions();

q0 = [0;0;.6;-1.2;.6+pi/2];
x0 = [q0;0*q0];
x0(2) = 1.0;
% x0=0.1*randn(r.getNumStates,1);
% x0(3) = 1.5;

tf = 0.5;

% Forward simulate dynamics with visulazation, then playback at realtime
true_traj = r.simulate([0 tf],x0);
true_traj = PPTrajectory(foh(true_traj.getBreaks(), true_traj.eval(true_traj.getBreaks())));
true_traj = true_traj.setOutputFrame(m.getStateFrame());
v.playback(true_traj,struct('slider',true));

N = 10:10:50;
rms_pos1 = zeros(length(N),1);
rms_pos2 = zeros(length(N),1);

for k = 1:length(N)
    N(k)
    dt = tf/(N(k)-1);
    ts = 0:dt:tf;
    
    r1 = TimeSteppingRigidBodyManipulator(m,dt,options);
    traj1{k} = r1.simulate([0 tf],x0);
    traj1{k} = PPTrajectory(foh(ts,traj1{k}.eval(traj1{k}.getBreaks)));
    traj1{k} = traj1{k}.setOutputFrame(m.getStateFrame());
    v.playback(traj1{k});
    
    opt = VariationalTrajectoryOptimization(m,N(k),tf,options);
    opt = opt.setSolver('snopt');
    opt = opt.setSolverOptions('snopt','MinorFeasibilityTolerance',1e-4);
    opt = opt.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-4);
    opt = opt.setSolverOptions('snopt','MinorOptimalityTolerance',1e-4);
    opt = opt.setSolverOptions('snopt','MajorOptimalityTolerance',1e-4);
    opt = opt.addPositionConstraint(ConstantConstraint(q0),1);
    opt = opt.addVelocityConstraint(ConstantConstraint(x0(6:end)),1);
    traj_init.x = true_traj;
    traj2{k} = opt.solveTraj(tf,traj_init);
    traj2{k} = traj2{k}.setOutputFrame(m.getStateFrame());
    v.playback(traj2{k});
    
    x1 = traj1{k}.eval(ts);
    x2 = traj2{k}.eval(ts);
    q_true = true_traj.eval(ts);
    rms_pos1(k) = rms(vec(x1(1:nq,:)-q_true(1:nq,:)));
    rms_pos2(k) = rms(vec(x2(1:nq,:)-q_true(1:nq,:)));
    
end

figure();
plot(rms_pos1);
hold on
plot(rms_pos2);
