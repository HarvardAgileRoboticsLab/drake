/// @file
///
/// kuka_plan_runner is designed to wait for LCM messages contraining
/// a robot_plan_t message, and then execute the plan on an iiwa arm
/// (also communicating via LCM using the
/// lcmt_iiwa_command/lcmt_iiwa_status messages).
///
/// When a plan is received, it will immediately begin executing that
/// plan on the arm (replacing any plan in progress).

#include <lcm/lcm-cpp.hpp>

#include "robotlocomotion/robot_plan_t.hpp" // This file is not found, to be figured out with other TRIers

#include "drake/common/drake_assert.h"
#include "drake/common/drake_path.h"
#include "drake/common/trajectories/piecewise_polynomial.h"
#include "drake/common/trajectories/piecewise_polynomial_trajectory.h"
#include "drake/examples/kuka_iiwa_arm/iiwa_common.h"
#include "drake/multibody/rigid_body_tree.h"

#include "drake/lcmt_iiwa_command.hpp"
#include "drake/lcmt_iiwa_status.hpp"
#include "drake/lcmt_polynomial.hpp" // temporarily abused lcm channel
#include "drake/lcmt_whole_body_data.hpp" // temporarily abused lcm channel

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
const char* const kLcmCommandChannel = "IIWA_COMMAND";
const char* const kLcmPlanChannel = "COMMITTED_ROBOT_PLAN";
const char* const kLcmParamChannel = "IIWA_PARAM";
const char* const kLcmParam2Channel = "IIWA_PARAM2";

const int kNumJoints = 7;

typedef PiecewisePolynomial<double> PPType;
typedef PPType::PolynomialType PPPoly;
typedef PPType::PolynomialMatrix PPMatrix;

class RobotPlanRunner {
 public:
  /// tree is aliased
  explicit RobotPlanRunner(const RigidBodyTree<double>& tree)
      : tree_(tree), plan_number_(0) {
    VerifyIiwaTree(tree);
    lcm_.subscribe(kLcmStatusChannel,
                    &RobotPlanRunner::HandleStatus, this);
    lcm_.subscribe(kLcmPlanChannel,
                    &RobotPlanRunner::HandlePlan, this);
  }

  void Run() {
    int cur_plan_number = plan_number_;
    int64_t cur_time_us = -1;
    int64_t start_time_us = -1;
    bool jointState_initialized = false;

    // Initialize the timestamp to an invalid number so we can detect
    // the first message.
    iiwa_status_.utime = cur_time_us;

    lcmt_iiwa_command iiwa_command;
    iiwa_command.num_joints = kNumJoints;
    iiwa_command.joint_position.resize(kNumJoints, 0.);
    iiwa_command.num_torques = kNumJoints;
    iiwa_command.joint_torque.resize(kNumJoints, 0.);

    lcmt_polynomial iiwa_param;
    iiwa_param.num_coefficients = kNumJoints;
    iiwa_param.coefficients.resize(kNumJoints, 0.);

    lcmt_whole_body_data iiwa_param2;
    iiwa_param2.num_positions = kNumJoints;
    iiwa_param2.q_des.resize(kNumJoints, 0.);
    iiwa_param2.num_constrained_dofs = 0;
    iiwa_param2.constrained_dofs.resize(0, 0);

    const int kNumDof = 7;
    Eigen::VectorXd joint_position_desired_previous(kNumDof); // 7DOF joint position at previous time sample
    Eigen::VectorXd joint_velocity_desired_previous(kNumDof); // 7DOF joint position at previous time sample

    while (true) {
      // Call lcm handle until at least one message is processed
      while (0 == lcm_.handleTimeout(10)) { }

      DRAKE_ASSERT(iiwa_status_.utime != -1);
      cur_time_us = iiwa_status_.utime;

      if (plan_) {
        if (plan_number_ != cur_plan_number) {
          std::cout << "Starting new plan." << std::endl;
          start_time_us = cur_time_us;
          cur_plan_number = plan_number_;
        }

        const double cur_traj_time_s =
            static_cast<double>(cur_time_us - start_time_us) / 1e6;
        const auto desired_next = plan_->value(cur_traj_time_s);

        Eigen::VectorXd jointPosPreviousState(kNumDof);
        if (!jointState_initialized){
          jointPosPreviousState << desired_next;
          jointState_initialized = true;
        }

        iiwa_command.utime = iiwa_status_.utime;
        iiwa_param.timestamp = iiwa_status_.utime;
        iiwa_param2.timestamp = iiwa_status_.utime;

        const double cur_time_inc_s = 0.005;
            
        // Inverse Dynamics
        // Set desired joint position and velocity
        Eigen::VectorXd joint_position_desired(kNumDof); // 7DOF joint position
        joint_position_desired.setZero();
        Eigen::VectorXd joint_velocity_desired(kNumDof); // 7DOF joint velocity
        joint_velocity_desired.setZero();
        Eigen::VectorXd joint_accel_desired(kNumDof); // 7DOF joint velocity
        joint_accel_desired.setZero();
        Eigen::VectorXd torque_ref(kNumDof);
        torque_ref.setZero();

        // Computing inverse dynamics torque command
        joint_position_desired << desired_next;
        KinematicsCache<double> cache = tree_.doKinematics(joint_position_desired, joint_velocity_desired);
        const RigidBodyTree<double>::BodyToWrenchMap no_external_wrenches;

        for (int joint = 0; joint < kNumJoints; joint++){
            joint_velocity_desired(joint) = (desired_next(joint) - joint_position_desired_previous(joint))/0.005;
            joint_accel_desired(joint) = (joint_velocity_desired(joint) - joint_velocity_desired_previous(joint))/0.005;
        }
          joint_position_desired_previous = desired_next;
          joint_velocity_desired_previous = joint_velocity_desired;

        // Use PD controller as the joint acceleration
        Eigen::VectorXd vd(kNumDof);
        for (int joint = 0; joint < kNumJoints; joint++) {
          vd(joint) = joint_accel_desired(joint);// velocity derivative
        }
        // note that, the gravity term in the inverse dynamics is set to zero.
        Eigen::VectorXd torque_command = torque_ref + tree_.inverseDynamics(cache, no_external_wrenches, vd, false);

        for (int joint = 0; joint < kNumJoints; joint++) {
            iiwa_param.coefficients[joint] = torque_command(joint);
        }

        // pure PD position control
        Eigen::VectorXd position_ctrl_torque_command(kNumDof);
        Eigen::VectorXd Kp_pos_ctrl(kNumDof); // 7 joints
        Kp_pos_ctrl << 225, 289, 144, 49, 324, 49, 49;
        Eigen::VectorXd Kd_pos_ctrl(kNumDof); // 7 joints
        Kd_pos_ctrl << 30, 34, 24, 14, 36, 14, 14;
        // (TODOs) Add integral control (anti-windup)
        for (int joint = 0; joint < kNumJoints; joint++) {
          position_ctrl_torque_command(joint) = Kp_pos_ctrl(joint)*(joint_position_desired(joint)
               - iiwa_status_.joint_position_measured[joint]) + Kd_pos_ctrl(joint)*(joint_velocity_desired(joint)
               - iiwa_status_.joint_velocity_estimated[joint]);
        }
        //Combination of ID torque control and IK position control
        torque_command += position_ctrl_torque_command;

        // -------->(For Safety) Set up iiwa position command<----------
        // Use the joint velocity estimation
        const double max_joint_velocity_estimated_term = 0.3; // [value to be optimized]
        for (int joint = 0; joint < kNumJoints; joint++) {
          double joint_velocity_estimated_term = iiwa_status_.joint_velocity_estimated[joint]*cur_time_inc_s;
          joint_velocity_estimated_term = std::max(-max_joint_velocity_estimated_term,
                                 std::min(max_joint_velocity_estimated_term, joint_velocity_estimated_term));
          iiwa_command.joint_position[joint] = iiwa_status_.joint_position_measured[joint] + joint_velocity_estimated_term;
        }

        // -------->Set up iiwa torque command<-------------
        for (int joint = 0; joint < kNumJoints; joint++) {
          iiwa_command.joint_torque[joint] = torque_command(joint);
          iiwa_command.joint_torque[joint] = std::max(-150.0,
                                                     std::min(150.0, iiwa_command.joint_torque[joint]));
          iiwa_param2.q_des[joint] = position_ctrl_torque_command(joint);
        }

        lcm_.publish(kLcmCommandChannel, &iiwa_command);
        lcm_.publish(kLcmParamChannel, &iiwa_param);
        lcm_.publish(kLcmParam2Channel, &iiwa_param2);
      }
    }
  }

 private:
  void HandleStatus(const lcm::ReceiveBuffer* rbuf, const std::string& chan,
                    const lcmt_iiwa_status* status) {
    iiwa_status_ = *status;
  }

  void HandlePlan(const lcm::ReceiveBuffer* rbuf, const std::string& chan,
                  const robotlocomotion::robot_plan_t* plan) {
    std::cout << "New plan received." << std::endl;
    Eigen::MatrixXd traj_mat(kNumJoints, plan->num_states);
    traj_mat.fill(0);

    std::map<std::string, int> name_to_idx =
        tree_.computePositionNameToIndexMap();
    for (int i = 0; i < plan->num_states; ++i) {
      const auto& state = plan->plan[i];
      for (int j = 0; j < state.num_joints; ++j) {
        if (name_to_idx.count(state.joint_name[j]) == 0) {
          continue;
        }
        traj_mat(name_to_idx[state.joint_name[j]], i) =
            state.joint_position[j];
      }
    }

    std::cout << traj_mat << std::endl;

    std::vector<double> input_time;
    for (int k = 0; k < static_cast<int>(plan->plan.size()); ++k) {
      input_time.push_back(plan->plan[k].utime / 1e6);
    }
    plan_.reset(new PiecewisePolynomialTrajectory(traj_mat, input_time));
    ++plan_number_;
  }

  lcm::LCM lcm_;
  const RigidBodyTree<double>& tree_;
  int plan_number_{};
  std::unique_ptr<PiecewisePolynomialTrajectory> plan_;
  lcmt_iiwa_status iiwa_status_;
};

int do_main(int argc, const char* argv[]) {
  RigidBodyTree<double> tree(
      drake::GetDrakePath() + "/examples/kuka_iiwa_arm/urdf/iiwa14.urdf",
      drake::multibody::joints::kFixed);

  RobotPlanRunner runner(tree);
  runner.Run();
  return 0;
}

}  // namespace
}  // namespace kuka_iiwa_arm
}  // namespace examples
}  // namespace drake


int main(int argc, const char* argv[]) {
  return drake::examples::kuka_iiwa_arm::do_main(argc, argv);
}
