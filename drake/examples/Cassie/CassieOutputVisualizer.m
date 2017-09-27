
options.use_bullet = false;
options.floating = 'quat';
r = Cassie(options);
v = r.constructVisualizer;

lc = lcm.lcm.LCM.getSingleton();
lcmAggregator = drake.matlab.util.MessageMonitor(drake.lcmt_cassie_state,'utime');
lc.subscribe('CASSIE_STATE',lcmAggregator)

motor_pos_to_state = [8 10 12 14 22 9 11 13 15 23]';
joint_pos_to_state = [16 18 17 19]'; 

nq = r.getNumPositions();
nv = r.getNumVelocities();
x = zeros(r.getNumStates(),1);

t=0;
dt = 0.01;
while true
  % check if there is a message available
  msg = lcmAggregator.getNextMessage(5);
  if isempty(msg)
    continue
  end
  
  msg = drake.lcmt_cassie_state(msg);

  q = msg.imu_orientation;
  q(1) = msg.imu_orientation(4);
  q(2:4) = msg.imu_orientation(1:3);
  qr=rotmat2quat(rotx(pi)); % IMU to pelvis transformation

  q = quatProduct(q,qr);


  % x(nq+(1:3)) = x(nq+(1:3)) + dt*(msg.imu_linear_acceleration - [0;0;9.81]);
  % x(1:3) = x(1:3) +dt*x(nq+(1:3)); % arb pelvis position for now
  x(3) = 1.0;
  x(4:7) = q;
  x(motor_pos_to_state) = msg.motor_positions;
  x(joint_pos_to_state) = msg.joint_positions;
  % x(nq+(4:6)) = msg.imu_angular_velocity;
  % x(nq+6:15) = msg.motor_velocities;

  t = t+0.01;
  v.draw(t,x);
  
  pause(0.05)
end