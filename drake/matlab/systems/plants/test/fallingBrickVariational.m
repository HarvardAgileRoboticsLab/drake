function fallingBrickVariational

options.floating = 'rpy';
options.terrain = RigidBodyFlatTerrain();
options.active_collision_options.terrain_only = true;
options.use_bullet = false;
s = 'FallingBrickContactPoints.urdf';
%p = TimeSteppingRigidBodyManipulator(s,.05,options);
p = VariationalRigidBodyManipulator3(s,.05,options);

%x0 = [0;1;2;randn(3,1);zeros(6,1)];
x0 = [0;0;1;pi/2-.2;.2;0;2;2;0;zeros(3,1)];

v = p.constructVisualizer();
v.drawWrapper(0,x0);
xtraj = p.simulate([0 1],x0);
v.playback(xtraj,'slider',true);