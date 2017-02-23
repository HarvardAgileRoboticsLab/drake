function fallingBall

options.twoD = 'true';
options.floating = 'true';
%options.floating = 'quat';
options.terrain = RigidBodyFlatTerrain();

s = 'Ball.urdf';
%p = TimeSteppingRigidBodyManipulator(s,.05,options);
p = VariationalTimeSteppingRigidBodyManipulator(s,.2,options);
x0 = [0;1;0;0;0;5];
%x0 = [0;1;2;rpy2quat(randn(3,1));randn(6,1)];
x0 = p.resolveConstraints(x0);
v = p.constructVisualizer();
v.drawWrapper(0,x0);
xtraj = p.simulate([0 3],x0);
v.playback(xtraj);

end