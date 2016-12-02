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

#include "robotlocomotion/robot_plan_t.hpp"

#include "drake/common/drake_assert.h"
#include "drake/common/drake_path.h"
#include "drake/common/trajectories/piecewise_polynomial.h"
#include "drake/common/trajectories/piecewise_polynomial_trajectory.h"
#include "drake/examples/kuka_iiwa_arm/iiwa_common.h"
#include "drake/multibody/rigid_body_tree.h"

#include "drake/lcmt_iiwa_status.hpp"
#include "drake/lcmt_robot_controller_reference.hpp"

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
const char* const kLcmPlanChannel = "COMMITTED_ROBOT_PLAN";
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

    // Initialize the timestamp to an invalid number so we can detect the first message.
    iiwa_status_.utime = cur_time_us;

    lcmt_robot_controller_reference robot_controller_reference;
    robot_controller_reference.num_joints = kNumJoints;
    robot_controller_reference.joint_position_desired.resize(kNumJoints, 0.);
    robot_controller_reference.joint_velocity_desired.resize(kNumJoints, 0.);
    robot_controller_reference.joint_accel_desired.resize(kNumJoints, 0.);
    robot_controller_reference.u_nominal.resize(kNumJoints, 0.);

    while (true) {
      // Call lcm handle until at least one message is processed
      while (0 == lcm_.handleTimeout(10)) { }

      DRAKE_ASSERT(iiwa_status_.utime != -1);
      cur_time_us = iiwa_status_.utime;

      if (qtraj_) {
        if (plan_number_ != cur_plan_number) {
          std::cout << "Starting new plan." << std::endl;
          start_time_us = cur_time_us;
          cur_plan_number = plan_number_;
        }
        const double cur_traj_time_s = static_cast<double>(cur_time_us - start_time_us) / 1e6;

        const auto q_ref = qtraj_->value(cur_traj_time_s);
        const auto qd_ref = qdtraj_->value(cur_traj_time_s);
        const auto qdd_ref = qddtraj_->value(cur_traj_time_s);

        robot_controller_reference.utime = iiwa_status_.utime;

        for(int joint = 0; joint < kNumJoints; joint++){
          robot_controller_reference.joint_position_desired[joint] = q_ref(joint);
          robot_controller_reference.joint_velocity_desired[joint] = qd_ref(joint);
          robot_controller_reference.joint_accel_desired[joint] = qdd_ref(joint);
        }

        lcm_.publish(kLcmControlRefChannel, &robot_controller_reference);
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
    qtraj_.reset(new PiecewisePolynomialTrajectory(traj_mat, input_time));
    qdtraj_.reset(new PiecewisePolynomialTrajectory(qtraj_->getPP().derivative(1)));
    qddtraj_.reset(new PiecewisePolynomialTrajectory(qdtraj_->getPP().derivative(1)));
    ++plan_number_;
  }

  lcm::LCM lcm_;
  const RigidBodyTree<double>& tree_;
  int plan_number_{};
  std::unique_ptr<PiecewisePolynomialTrajectory> qtraj_;
  std::unique_ptr<PiecewisePolynomialTrajectory> qdtraj_;
  std::unique_ptr<PiecewisePolynomialTrajectory> qddtraj_;
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