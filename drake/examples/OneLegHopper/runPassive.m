warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
warning('off','Drake:RigidBodyManipulator:WeldedLinkInd');
warning('off','Drake:RigidBodyManipulator:UnsupportedJointLimits');
options.twoD = true;
options.view = 'right';
options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.ignore_self_collisions = true;
s = 'OneLegHopper.urdf';
dt = 0.005;
r = TimeSteppingRigidBodyManipulator(s,dt,options);
r = r.setStateFrame(OneLegHopperState(r));
r = r.setOutputFrame(OneLegHopperState(r));
v = r.constructVisualizer();

q0 = [0;0;.6;-1.2;.6+pi/2];
x0 = [q0;0*q0];
x0(2) = 1.0;
% x0=0.1*randn(r.getNumStates,1);
% x0(3) = 1.5;

% Forward simulate dynamics with visulazation, then playback at realtime
S=warning('off','Drake:DrakeSystem:UnsupportedSampleTime');
output_select(1).system=1;
output_select(1).output=1;
sys = mimoCascade(r,v,[],[],output_select);
warning(S);
traj = simulate(sys,[0 2],x0);
playback(v,traj,struct('slider',true));

