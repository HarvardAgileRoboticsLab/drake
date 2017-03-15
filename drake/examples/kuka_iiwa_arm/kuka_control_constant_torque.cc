/// @file
///
/// kuka_control is designed to compute the torque command based on 
/// desired joint position, velocity, and acceleration, and measured joint position and velocity.
/// Currently, the final torque command is composed of inverse dynamics torque command and joint position PD
/// controller command.
/// (TODO: Generalize this ID controller to more general types of feedback controllers)
/// Messages are sent via LCM channels.

#include <lcm/lcm-cpp.hpp>

#include "drake/common/drake_assert.h"
#include "drake/common/drake_path.h"
#include "drake/examples/kuka_iiwa_arm/iiwa_common.h"
#include "drake/multibody/parsers/urdf_parser.h"
#include "drake/multibody/rigid_body_tree.h"

#include "drake/lcmt_robot_controller_reference.hpp"
#include "drake/lcmt_iiwa_status.hpp"
#include "drake/lcmt_iiwa_command.hpp"

#include "drake/util/drakeGeometryUtil.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using drake::Vector1d;
using Eigen::Vector2d;
using Eigen::Vector3d;

namespace drake {
namespace examples {
namespace kuka_iiwa_arm {
namespace {

const char* const kLcmStatusChannel = "IIWA_STATUS";
const char* const kLcmControlRefChannel = "CONTROLLER_REFERENCE";
const char* const kLcmCommandChannel = "IIWA_COMMAND";
const int kNumJoints = 7;

class RobotController {
 public:
  /// tree is aliased
  explicit RobotController(const RigidBodyTree<double>& tree)
      : tree_(tree), controller_trigger_(false) {
    VerifyIiwaTree(tree);
    lcm_.subscribe(kLcmStatusChannel,
                    &RobotController::HandleStatus, this);
    lcm_.subscribe(kLcmControlRefChannel,
                    &RobotController::HandleControl, this);
  }

  void Run() {
    int64_t cur_time_us = -1;

    // Initialize the timestamp to an invalid number so we can detect the first message.
    iiwa_status_.utime = cur_time_us;
    robot_controller_reference_.utime = cur_time_us;

    lcmt_iiwa_command iiwa_command;
    iiwa_command.num_joints = kNumJoints;
    iiwa_command.joint_position.resize(kNumJoints, 0.);
    iiwa_command.num_torques = kNumJoints;
    iiwa_command.joint_torque.resize(kNumJoints, 0.);

    Eigen::VectorXd joint_position_desired(kNumJoints); 
    Eigen::VectorXd joint_velocity_desired(kNumJoints); 
    Eigen::VectorXd joint_accel_desired(kNumJoints); 

    bool half_servo_rate_flag_ = true; // make the iiwa command get published every two servo loops

    while (true) {
      // Call lcm handle until at least one message is processed
      while (0 == lcm_.handleTimeout(10)) { }
      DRAKE_ASSERT(iiwa_status_.utime != -1);
      DRAKE_ASSERT(robot_controller_reference_.utime != -1);
      cur_time_us = iiwa_status_.utime;

      if (controller_trigger_) {  
        const int kNumDof = 7;
        iiwa_command.utime = iiwa_status_.utime;

        // Kuka-Controller (Currently implement an inverse dynamics controller)
        // Set desired joint position, velocity and acceleration
        for (int joint = 0; joint < kNumDof; joint++){
          joint_position_desired(joint) = robot_controller_reference_.joint_position_desired[joint];
          joint_velocity_desired(joint) = robot_controller_reference_.joint_velocity_desired[joint];
          joint_accel_desired(joint) = robot_controller_reference_.joint_accel_desired[joint];
        }

        double *qptr = &iiwa_status_.joint_position_measured[0];
        Eigen::Map<Eigen::VectorXd> q(qptr, kNumDof);
        double *qdptr = &iiwa_status_.joint_velocity_estimated[0];
        Eigen::Map<Eigen::VectorXd> qd(qdptr, kNumDof);

        Eigen::VectorXd torque_command(kNumDof);
        torque_command.setZero();
        //torque_command(4) = 0;// to be tested
        //torque_command(0) = 3.3;
        //torque_command(1) = 5;
        //torque_command(2) = 1.65;
        //torque_command(3) = -1.39; //[0 --> positive position: 1.43 Nm, 0 --> negative position: -1.36Nm, 
                                     // negative position --> 0: 1.40 Nm, positive position --> 0: -1.39 Nm]

        // -------->(For Safety) Set up iiwa position command<----------
        for (int joint = 0; joint < kNumJoints; joint++) {
          iiwa_command.joint_position[joint] = joint_position_desired(joint);
        }

        // -------->Set up iiwa torque command<-------------
        for (int joint = 0; joint < kNumJoints; joint++) {
          iiwa_command.joint_torque[joint] = torque_command(joint);
          iiwa_command.joint_torque[joint] = std::max(-150.0, std::min(150.0, iiwa_command.joint_torque[joint]));
        }

        if (half_servo_rate_flag_){
          half_servo_rate_flag_ = false;
          lcm_.publish(kLcmCommandChannel, &iiwa_command);
        }else{
          half_servo_rate_flag_ = true;
        }
      }
    
    }
  }

 private:
  void HandleStatus(const lcm::ReceiveBuffer* rbuf, const std::string& chan,
                    const lcmt_iiwa_status* status) {
    iiwa_status_ = *status;
  }

  void HandleControl(const lcm::ReceiveBuffer* rbuf, const std::string& chan,
                    const lcmt_robot_controller_reference* input) {
    robot_controller_reference_ = *input;
    controller_trigger_ = true;
  }

  lcm::LCM lcm_;
  const RigidBodyTree<double>& tree_;
  bool controller_trigger_;// control runner wait for the first message from plan runner 
  lcmt_iiwa_status iiwa_status_;
  lcmt_robot_controller_reference robot_controller_reference_;
};

int DoMain(int argc, const char* argv[]) {

  auto tree = std::make_unique<RigidBodyTree<double>>();
  parsers::urdf::AddModelInstanceFromUrdfFileToWorld(
      GetDrakePath() + "/examples/kuka_iiwa_arm/urdf/iiwa14.urdf",
      multibody::joints::kFixed, tree.get());

  RobotController runner(*tree);
  runner.Run();
  return 0;
}

}  // namespace
}  // namespace kuka_iiwa_arm
}  // namespace examples
}  // namespace drake

int main(int argc, const char* argv[]) {
  return drake::examples::kuka_iiwa_arm::DoMain(argc, argv);
}