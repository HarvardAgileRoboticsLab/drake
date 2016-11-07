#include <iostream>
#include <time.h>
#include <math.h>

#include "drake/examples/kuka_iiwa_arm/kuka_message_handler.h"

#include "drake/systems/plants/KinematicsCache.h"

#include "drake/common/drake_assert.h"
#include "drake/common/drake_path.h"
#include "drake/common/polynomial.h"
#include "drake/examples/kuka_iiwa_arm/iiwa_status.h"
#include "drake/systems/plants/IKoptions.h"
#include "drake/systems/plants/RigidBodyIK.h"
#include "drake/systems/plants/RigidBodyTree.h"
#include "drake/systems/plants/constraint/RigidBodyConstraint.h"
#include "drake/systems/plants/joints/floating_base_types.h"
#include "drake/systems/trajectories/PiecewisePolynomial.h"
#include "drake/systems/vector.h"

namespace drake {
namespace examples {
namespace kuka_iiwa_arm {
namespace{


const int kNumDOF = 7;

int main(int argc, const char* argv[]){
  // allocate a pointer to an LCM object and create the message handler
  KukaMessageHandler handler;

  bool debug = true;
  if (!debug){// for debugging
    // wait for the first status to come in
    bool success = handler.waitForFirstMessage();

    if (!success){
      printf("Timed out waiting for status\n");
      return 1;
    }
  }

  // load the model for the kuka arm
  RigidBodyTree<double> tree(
    drake::GetDrakePath() + "/examples/kuka_iiwa_arm/urdf/iiwa14.urdf",
    drake::systems::plants::joints::kFixed);

  // initialize inverse dynamics inputs
  Eigen::VectorXd position;
  if (debug){
    position = tree.getZeroConfiguration();
  }else{
    position = handler.getPosition();
  }
  Eigen::VectorXd accelerations = Eigen::VectorXd::Constant(kNumDOF,1,0.0);
  // initialize inverse dynamcis output
  
  // setup kinematics cache 
  KinematicsCache<double> kinCache(tree.bodies);
  kinCache.initialize(position);
  tree.doKinematics(kinCache);
  RigidBodyTree<double>::BodyToWrenchMap no_external_wrenches;

// setup the inverse kinematics
  const int kNumTimesteps = 20;
  const double dt = 0.1;
  double t[kNumTimesteps];
  for (int i=0; i < kNumTimesteps; i++){
    t[i] = i*dt;
  }
  Eigen::MatrixXd q0(tree.get_num_positions(), kNumTimesteps);
  Eigen::VectorXd qdot0 = Eigen::VectorXd::Zero(tree.get_num_velocities());
  Eigen::MatrixXd q_sol(tree.get_num_positions(), kNumTimesteps);
  Eigen::MatrixXd qdot_sol(tree.get_num_velocities(), kNumTimesteps);
  Eigen::MatrixXd qddot_sol(tree.get_num_positions(), kNumTimesteps);
  Eigen::MatrixXd torques(tree.get_num_positions(), kNumTimesteps);
  for (int i = 0; i < kNumTimesteps; i++) {
    q0.col(i) = position;
  }
  IKoptions ik_options(&tree);
  ik_options.setMajorIterationsLimit(50000);

  // construct the constraints, starting at the current configuration
  std::vector<RigidBodyConstraint*> constraint_array;
  
  Eigen::VectorXd joint_eps = Eigen::VectorXd::Constant(7, 0.1);
  Eigen::Vector3d world_eps = Eigen::Vector3d::Constant(0.1);

  // the current position of the robot
  PostureConstraint pc0(&tree, Eigen::Vector2d(0, 0.5));
  Eigen::VectorXi joint_idx(7);
  joint_idx << 0, 1, 2, 3, 4, 5, 6;
  Eigen::VectorXd joint_lb = position - joint_eps;
  Eigen::VectorXd joint_ub = position + joint_eps;
  pc0.setJointLimits(joint_idx, joint_lb, joint_ub);
  constraint_array.push_back(&pc0);

  // the desired location
  Eigen::Vector3d pos_end;
  if (argc <= 3){
    pos_end << 0.6, 0, 0.325;
  }else{
    pos_end << atof(argv[1]), atof(argv[2]), atof(argv[3]); 
  }
  printf("Navigating the end effector to %f, %f, %f\n",pos_end[0],pos_end[1],pos_end[2]);
  Eigen::Vector3d pos_lb = pos_end - world_eps;
  Eigen::Vector3d pos_ub = pos_end + world_eps;
  WorldPositionConstraint wpc(&tree, tree.FindBodyIndex("iiwa_link_ee"),
                              Eigen::Vector3d::Zero(), pos_lb, pos_ub, Eigen::Vector2d(t[kNumTimesteps-1], t[kNumTimesteps-1]));
  constraint_array.push_back(&wpc);

  int info = 0;
  std::vector<std::string> infeasible_constraint;
  inverseKinTraj(&tree, kNumTimesteps, t, qdot0, q0, q0, constraint_array.size(),
                constraint_array.data(), ik_options, &q_sol, &qdot_sol, &qddot_sol,
                &info, &infeasible_constraint);

  printf("INFO = %d ", info);
  printf("\n");

  if (info >= 10) {
    std::cerr << "Solution failed, not sending." << std::endl;
    return 1;
  }


  // use the inverse dynamics to compute the torques along the trajectory
  Eigen::VectorXd pos, vel, acc;
  for (int i=0; i<kNumTimesteps; i++){
    pos = q_sol.col(i);
    vel = qdot_sol.col(i);
    acc = qddot_sol.col(i);

    kinCache = tree.doKinematics(pos,vel);
    
    torques.col(i) = tree.inverseDynamics(kinCache,  
                      no_external_wrenches, acc, false);

    std::cout << torques.col(i).transpose() << std::endl;
  }

  // run the sequence
  double t_init, t_now;
  if (debug){
    t_init = time(NULL);
  }else{
    t_init = handler.getTime();
  }
  t_now = t_init;
  // control loop
  Eigen::VectorXd pos_measured(tree.get_num_positions());
  Eigen::VectorXd vel_measured(tree.get_num_positions());
  Eigen::VectorXd input(tree.get_num_positions());
  Eigen::VectorXd err(tree.get_num_positions());
  Eigen::VectorXd G(tree.get_num_positions());
  int t_idx=0;
  while(!handler.hasTimedOut() || (debug && t_now-t_init < t[kNumTimesteps-1])){
    if (handler.handle() || debug){// if a new state has come in 
      if(debug){
        t_now = time(NULL);
      }else{
        t_now = handler.getTime();
      }
      
      t_idx = floor((t_now-t_init)/dt);

      // torque from plan
      input = torques.col(t_idx);

      // implement proportional control
      if (debug){
        pos_measured = q_sol.col(t_idx);
        vel_measured = qdot_sol.col(t_idx);
      }else{
        pos_measured = handler.getPosition();
        vel_measured = handler.getVelocity();
      }

      err = pos_measured - q_sol.col(t_idx);

      // compute the gravity compensation
      kinCache = tree.doKinematics(pos_measured, vel_measured);
      G = tree.inverseDynamics(kinCache, no_external_wrenches,
                               accelerations, false);

      input = input - 0.1*err + G; // TODO: make this tvlqr
      
      // publish the state
      if (debug){
        printf("publishing torque\n");
        std::cout << input.transpose() << std::endl;
      }else{
        handler.publishTorques(-1*torques);
      }
    }
  }

  if(!handler.hasTimedOut()){
    printf("Timed out waiting for status\n");
    return 1;
  }

  return 0;
}

} // namespace
} // kuka_iiwa_arm
} // examples
} // drake

int main(int argc, const char* argv[]) {
  return drake::examples::kuka_iiwa_arm::main(argc, argv);
}