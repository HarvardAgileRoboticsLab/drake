function fallingBall

%options.twoD = 'true';
options.enable_fastqp = false;
options.floating = true;
%options.floating = 'quat';
options.terrain = RigidBodyFlatTerrain();

s = 'ball.urdf';

%x0 = [0;1;0;0;0;10];
x0 = [0;0;1;0;0;0;0;2;0;0;0;0];

p1 = TimeSteppingRigidBodyManipulator(s,.05,options);
x01 = p1.resolveConstraints(x0);
v1 = p1.constructVisualizer();
v1.drawWrapper(0,x0);
tic
xtraj1 = p1.simulate([0 3],x01);
toc
v1.playback(xtraj1);

p2 = VariationalRigidBodyManipulator(s,.1,options);
x02 = p2.resolveConstraints(x0);
tic
xtraj2 = p2.simulate([0 3],x02);
toc
v2 = p2.constructVisualizer();
v2.playback(xtraj2);

end