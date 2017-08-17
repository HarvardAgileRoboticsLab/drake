function fallingBrickLCP

options.floating = true;
options.terrain = RigidBodyFlatTerrain();
% options.ignore_self_collisions = true;
% options.use_bullet = false;
s = 'FallingBrickContactPoints.urdf';
% s = 'FallingBrickBetterCollisionGeometry.urdf';
p = TimeSteppingRigidBodyManipulator(s,.01,options);
% p = p.addRobotFromURDF(s,[],[],options);
x0 = [0;0;1;zeros(3,1);10;zeros(5,1)];
% x0 = [0;1;2;rpy2quat(randn(3,1));2;1;2;rpy2quat(randn(3,1));randn(12,1)];
x0 = p.resolveConstraints(x0);

if 0 
  v = p.constructVisualizer();
  sys = cascade(p,v);
  sys.simulate([0 8],x0);
  return;
end

v = p.constructVisualizer();
v.drawWrapper(0,x0);

% Forward simulate dynamics with visulazation, then playback at realtime
S=warning('off','Drake:DrakeSystem:UnsupportedSampleTime');
output_select(1).system=1;
output_select(1).output=1;
sys = mimoCascade(p,v,[],[],output_select);
warning(S);
traj = simulate(sys,[0 4],x0);
playback(v,traj,struct('slider',true));

