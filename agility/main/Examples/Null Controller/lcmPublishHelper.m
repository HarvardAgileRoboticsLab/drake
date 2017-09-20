function lcmPublishHelper(cassie_output)

  global lc;

  if isempty(lc)
    lc = lcm.lcm.LCM.getSingleton();  
  end

  state = drake.lcmt_cassie_state();

  state.utime = 0;
  state.status = cassie_output.status;
  state.charge_state = cassie_output.stateOfCharge;
  state.radio = cassie_output.radio;
  
  state.joint_positions = cassie_output.jointPosition;
  state.joint_velocities = cassie_output.jointVelocity;

  state.motor_positions = cassie_output.motorPosition;
  state.motor_velocities = cassie_output.motorVelocity;

  state.imu_angular_velocity = cassie_output.vectorNavAngularVelocity;
  state.imu_linear_acceleration = cassie_output.vectorNavLinearAcceleration;
  state.imu_orientation = cassie_output.vectorNavOrientation;
  state.imu_magnetic_field = cassie_output.vectorNavMagneticField;

  lc.publish('CASSIE_STATE', state);

end

