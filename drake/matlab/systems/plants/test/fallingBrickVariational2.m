function fallingBrickVariational2

options.floating = 'rpy';
options.terrain = RigidBodyFlatTerrain();
options.use_bullet = false;
s = 'FallingBrickContactPoints.urdf';
%p = TimeSteppingRigidBodyManipulator(s,.01,options);
p = VariationalRigidBodyManipulator2(s,.05,options);

%x0 = [0;1;2;randn(3,1);zeros(6,1)];
x0 = [0;0;1.5;pi/2-.2;.2;0;5;5;0;zeros(3,1)];

v = p.constructVisualizer();
v.drawWrapper(0,x0);
xtraj = p.simulate([0 2],x0);
v.playback(xtraj);