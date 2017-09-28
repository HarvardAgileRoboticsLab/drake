function testPrismaticSpring 

r = RigidBodyManipulator('PrismaticSpring.urdf');
v = r.constructVisualizer(); 

% note: these have to match the urdf
m = 1;

k = 10;
b = 1;
g = 9.81;
rest_angle = 0;

xtraj = r.simulate([0, 10], [0;0]);
v.playback(xtraj)

for i=1:100
  x = randn(2,1);
  xdot = dynamics(r,0,x,[]);
  theta_ddot_desired = (-b*x(2) -m*g + k*(rest_angle - x(1)))/m;
  valuecheck(xdot(2),theta_ddot_desired);
end

prismatic_spring = r.force{1};
q = getRandomConfiguration(r);
qd = rand(r.getNumVelocities());

geval_options.grad_method = {'user', 'taylorvar'};
[~, ~] = geval(1, @(q, qd) prismatic_spring.computeSpatialForce(r, q, qd), q, qd, geval_options);

end
