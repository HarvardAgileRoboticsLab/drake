function fallingBrickLCP

options.floating = 'quat';
options.terrain = RigidBodyFlatTerrain();
% options.ignore_self_collisions = true;
% options.use_bullet = false;
s = 'FallingBrickContactPoints.urdf';
% s = 'FallingBrickBetterCollisionGeometry.urdf';
p = TimeSteppingRigidBodyManipulator_energy_based_method(s,.01,options);
p = p.addRobotFromURDF(s,[],[],options);
% x0 = [0;1;2;rpy2quat(randn(3,1));randn(6,1)];

%sampled data set 1
fix1 = [0.4400; 0.1017; 2.7873];
fix2 = [-1.1667; -1.8543; -1.1407];
fix3 = [-0.2185    0.5413    0.3893    0.7512    1.7783    1.2231   -1.2833   -2.3290    0.9019   -1.8356    0.0668    0.0355]';

% %sampled data set 2
% fix1 = [   -0.0375; -1.8963; -2.1280];
% fix2 = [-1.1769; -0.9905; -1.1730];
% fix3 = [-1.7254 0.2882 -1.5942 0.1102 0.7871 -0.0022 0.0931 -0.3782 -1.4827 -0.0438 0.9608 1.7382]';

% x0 = [0;1;2;rpy2quat(randn(3,1));2;1;2;rpy2quat(randn(3,1));randn(12,1)];
x0 = [0;1;2;rpy2quat(fix1);2;1;2;rpy2quat(fix2);fix3];
x0 = p.resolveConstraints(x0);

if 0 
  v = p.constructVisualizer();
  sys = cascade(p,v);
  sys.simulate([0 8],x0);
  return;
end

v = p.constructVisualizer();
v.drawWrapper(0,x0);
xtraj = p.simulate([0 8],x0);
v.playback(xtraj);