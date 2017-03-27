urdf = [getDrakePath(), '/matlab/systems/plants/test/FallingBrick.urdf'];
options.floating = true;
options.terrain = RigidBodyFlatTerrain;

w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
r = RigidBodyManipulator(urdf,options);
r = r.addGeometryToBody('world',RigidBodySphere(1),'groupA');
r = compile(r);
warning(w);

x0 = Point(r.getStateFrame());
x0.base_z = 2;


r.constructVisualizer

