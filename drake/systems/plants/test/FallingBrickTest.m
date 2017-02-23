p = TimeSteppingRigidBodyManipulator('FallingBrick.urdf', .01);

x0 = [0 0 1 0 0 0 0 0 0 0 0 0]';
x0 = p.resolveConstraints(x0);
v = p.constructVisualizer();
v.drawWrapper(0,x0);
xtraj = p.simulate([0 1],x0);
v.playback(xtraj);