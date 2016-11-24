#include <iostream>
#include <time.h>

#include <lcm/lcm-cpp.hpp>

#include "drake/lcmt_iiwa_command.hpp"
#include "drake/lcmt_iiwa_status.hpp"

#include "drake/examples/kuka_iiwa_arm/iiwa_status.h"
#include "drake/examples/kuka_iiwa_arm/kuka_message_handler.h"
#include "drake/systems/plants/RigidBodyTree.h"
#include "drake/common/drake_path.h"
#include "drake/systems/plants/KinematicsCache.h"


namespace drake {
namespace examples {
namespace kuka_iiwa_arm {
namespace{


const int kNumDOF = 7;

int main(int argc, const char* argv[]){
  // allocate a pointer to an LCM object and create the message handler
  std::shared_ptr<lcm::LCM> lcm = std::make_shared<lcm::LCM>();
  KukaMessageHandler handler(lcm);

  // load the model for the kuka arm
  RigidBodyTree tree(
    drake::GetDrakePath() + "/examples/kuka_iiwa_arm/urdf/iiwa14.urdf",
    drake::systems::plants::joints::kFixed);
  
  Eigen::VectorXd position = Eigen::VectorXd::Constant(kNumDOF,1,0.0);
  Eigen::VectorXd accelerations = Eigen::VectorXd::Constant(kNumDOF,1,0.0);
  Eigen::VectorXd torques;
  KinematicsCache<double> kinCache(tree.bodies);
  kinCache.initialize(position);
  tree.doKinematics(kinCache);
  RigidBodyTree::BodyToWrenchMap<double> no_external_wrenches;

  // wait for the first status to come in
  bool success = handler.waitForFirstMessage();

  if (!success){
    printf("Timed out waiting for status\n");
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