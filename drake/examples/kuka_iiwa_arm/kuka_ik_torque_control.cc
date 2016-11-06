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
  RigidBodyTree tree(
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
  Eigen::VectorXd torques;
  // setup kinematics cache 
  KinematicsCache<double> kinCache(tree.bodies);
  kinCache.initialize(position);
  tree.doKinematics(kinCache);
  RigidBodyTree::BodyToWrenchMap<double> no_external_wrenches;

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
  if (argc <= 2){
    pos_end << 0.6, 0, 0.325;
  }else{
    pos_end << atof(argv[0]), atof(argv[1]), atof(argv[2]); 
  }
  Eigen::Vector3d pos_lb = pos_end - world_eps;
  Eigen::Vector3d pos_ub = pos_end + world_eps;
  WorldPositionConstraint wpc(&tree, tree.FindBodyIndex("iiwa_link_ee"),
                              Eigen::Vector3d::Zero(), pos_lb, pos_ub, Eigen::Vector2d(3.5, 4));
  constraint_array.push_back(&wpc);

  // setup the inverse kinematics
  const int kNumTimesteps = 20;
  Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(20,0,4);
  Eigen::MatrixXd q0(tree.get_num_positions(), kNumTimesteps);
  for (int i = 0; i < kNumTimesteps; i++) {
    q0.col(i) = position;
  }
  IKoptions ik_options(&tree);

  IKResults result = inverseKinTrajSimple(&tree, t, q0, q0, constraint_array, ik_options);
  
  bool info_good = true;
  for (int i = 0; i < kNumTimesteps; ++i) {
    printf("INFO[%d] = %-3d ", i, result.info[i]);
    switch (result.info[i]) {
      case 1:
        printf("Solution Found\n"); break;
      case 91:
        printf("Invalid Input\n"); break;
      case 13:
        printf("Infeasible Constraint\n"); break;
      default:// probably 100
        printf("Unknown Error\n"); break;
    }
    if (result.info[i] != 1) {
      info_good = false;
    }
  }
  printf("\n");


  // Debugging
  for (int i=0; i<result.q_sol.size(); i++){
    std::cout << result.q_sol[i] << std::endl;
  }

  if (!info_good) {
    std::cerr << "Solution failed, not sending." << std::endl;
    return 1;
  }

  // control loop
  while(!handler.hasTimedOut()){
    if (handler.handle()){// if a new state has come in

      // get the current position
      position = handler.getPosition();
      // compute the inverse dynamics
      kinCache = tree.doKinematics(position);
      torques = tree.inverseDynamics(kinCache, 
                      no_external_wrenches,
                      accelerations,
                      false); 
      
      // publish the state
      handler.publishTorques(-1*torques);
    }
  }
  printf("Timed out waiting for status\n");

  return 0;
}

} // namespace
} // kuka_iiwa_arm
} // examples
} // drake

int main(int argc, const char* argv[]) {
  return drake::examples::kuka_iiwa_arm::main(argc, argv);
}