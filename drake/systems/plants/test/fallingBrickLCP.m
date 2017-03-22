function fallingBrickLCP

%options.floating = 'quat';
%options.twoD = true;
options.floating = true;
options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.use_bullet = false;

s = 'FallingBrickContactPoints2.urdf';
%s = 'FallingBrick.urdf';
%s = 'FallingBrickBetterCollisionGeometry.urdf';

%p = TimeSteppingRigidBodyManipulator(s,.05,options);
p = VariationalTimeSteppingRigidBodyManipulator(s,.1,options);

%p = p.addRobotFromURDF(s,[],[],options);
%x0 = [0;1;2;rpy2quat(randn(3,1));randn(6,1)];
%x0 = [0;1;2;rpy2quat(randn(3,1));2;1;2;rpy2quat(randn(3,1));randn(12,1)];
x0 = [0;0;1;.1;.1;0;0;0;0;0;0;0];
%x0 = [0;2;0;0;0;0];
x0 = p.resolveConstraints(x0);

v = p.constructVisualizer();
v.drawWrapper(0,x0);

tic
xtraj = p.simulate([0 2],x0);
toc
v.playback(xtraj);