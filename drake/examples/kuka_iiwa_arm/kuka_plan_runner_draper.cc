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

#include "drake/lcmt_iiwa_command.hpp"
#include "drake/lcmt_iiwa_status.hpp"

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
const char* const kCancelPlanRunning = "CANCEL_PLAN";
const char* const kLcmGravityChannel = "IIWA_GRAVITY_TORQUE";

const int kNumJoints = 7;

typedef PiecewisePolynomial<double> PPType;
typedef PPType::PolynomialType PPPoly;
typedef PPType::PolynomialMatrix PPMatrix;

class RobotPlanRunner {
 public:
   bool run_ = false;
  /// tree is aliased
  explicit RobotPlanRunner(const RigidBodyTree<double>& tree)
      : tree_(tree), plan_number_(0) {
    VerifyIiwaTree(tree);
    lcm_.subscribe(kLcmStatusChannel,
                    &RobotPlanRunner::HandleStatus, this);
    lcm_.subscribe(kLcmPlanChannel,
                    &RobotPlanRunner::HandlePlan, this);
    lcm_.subscribe(kCancelPlanRunning,
                    &RobotPlanRunner::HandleCancelPlan, this);
  }

  void Run() {
    int cur_plan_number = plan_number_;
    int64_t cur_time_us = -1;
    int64_t start_time_us = -1;

    // Initialize the timestamp to an invalid number so we can detect
    // the first message.
    iiwa_status_.utime = cur_time_us;

    lcmt_iiwa_command iiwa_command;
    iiwa_command.num_joints = kNumJoints;
    iiwa_command.joint_position.resize(kNumJoints, 0.);
    iiwa_command.num_torques = kNumJoints;
    iiwa_command.joint_torque.resize(kNumJoints, 0.);

    while (true) {
      // Call lcm handle until at least one message is processed
      while (0 == lcm_.handleTimeout(10)) { }



      // const int kNumDof = 7;
      //
      // double *qptr = &iiwa_status_.joint_position_measured[0];
      // Eigen::Map<Eigen::VectorXd> q(qptr, kNumDof);
      // double *qdptr = &iiwa_status_.joint_velocity_estimated[0];
      // Eigen::Map<Eigen::VectorXd> qd(qdptr, kNumDof);
      // Eigen::VectorXd z = Eigen::VectorXd::Zero(kNumDof);
      //
      // // Computing inverse dynamics torque command
      // KinematicsCache<double> cache = tree_.doKinematics(q, qd);
      // const RigidBodyTree<double>::BodyToWrenchMap no_external_wrenches;
      // Eigen::VectorXd torque_command = tree_.inverseDynamics(cache, no_external_wrenches, joint_accel_desired, false);
      // // gravity compensation without gripper (to cancel out the low-level kuka controller)
      // Eigen::VectorXd z = Eigen::VectorXd::Zero(kNumDof);
      // Eigen::VectorXd gravity_torque = gravity_comp_no_gripper(cache, z, false, tree_);
      //




      if (!run_){
        // std::cout << "I am quitting!  \n" << run_ << std::endl;
        iiwa_status_.utime = iiwa_status_.utime;

        lcmt_iiwa_command iiwa_command;
        iiwa_command.num_joints = kNumJoints;
        iiwa_command.joint_position.resize(kNumJoints, 0.);
        iiwa_command.num_torques = kNumJoints;
        iiwa_command.joint_torque.resize(kNumJoints, 0.);
      }

      DRAKE_ASSERT(iiwa_status_.utime != -1);
      cur_time_us = iiwa_status_.utime;

      if (qtraj_) {
        if (plan_number_ != cur_plan_number) {
          std::cout << "Starting new plan." << std::endl;
          start_time_us = cur_time_us;
          cur_plan_number = plan_number_;
        }

        const int kNumDof = 7;
        const double cur_traj_time_s = static_cast<double>(cur_time_us - start_time_us) / 1e6;

        const auto q_ref = qtraj_->value(cur_traj_time_s);
        const auto qd_ref = qdtraj_->value(cur_traj_time_s);
        const auto qdd_ref = qddtraj_->value(cur_traj_time_s);

        iiwa_command.utime = iiwa_status_.utime;

        // Inverse Dynamics
        // Set desired joint position and velocity
        Eigen::VectorXd joint_position_desired(kNumDof);
        joint_position_desired << q_ref;
        Eigen::VectorXd joint_velocity_desired(kNumDof);
        joint_velocity_desired << qd_ref;
        Eigen::VectorXd joint_accel_desired(kNumDof);
        joint_accel_desired << qdd_ref;


        double *qptr = &iiwa_status_.joint_position_measured[0];
        Eigen::Map<Eigen::VectorXd> q(qptr, kNumDof);
        double *qdptr = &iiwa_status_.joint_velocity_estimated[0];
        Eigen::Map<Eigen::VectorXd> qd(qdptr, kNumDof);

        // Computing inverse dynamics torque command
        KinematicsCache<double> cache = tree_.doKinematics(q, qd);
        const RigidBodyTree<double>::BodyToWrenchMap no_external_wrenches;

        // note that, the gravity term in the inverse dynamics is set to zero.
        Eigen::VectorXd torque_command = tree_.inverseDynamics(cache, no_external_wrenches, joint_accel_desired, false);

        Eigen::VectorXd z = Eigen::VectorXd::Zero(kNumDof);
        KinematicsCache<double> cache0 = tree_.doKinematics(q, qd);
        // Eigen::VectorXd gravity_torque = tree_.inverseDynamics(cache0, no_external_wrenches, z, false);
        Eigen::VectorXd gravity_torque = gravity_comp_no_gripper(cache, z, false, tree_);

        torque_command -= gravity_torque;

        lcmt_iiwa_status grav_status;
        grav_status.utime = -1;
        grav_status.num_joints = kNumJoints;
        grav_status.joint_position_measured.resize(kNumJoints);
        grav_status.joint_position_commanded.resize(kNumJoints);
        grav_status.joint_position_ipo.resize(kNumJoints);
        grav_status.joint_velocity_estimated.resize(kNumJoints);
        grav_status.joint_torque_measured.resize(kNumJoints);
        grav_status.joint_torque_commanded.resize(kNumJoints);
        grav_status.joint_torque_external.resize(kNumJoints);
        grav_status.joint_acceleration_estimated.resize(kNumJoints);
        // grav_status.joint_torque_measured.resize(kNumJoints, 0.);
        // grav_status.utime = 0;
        for(int i = 0; i < kNumJoints; i++) {
          grav_status.joint_torque_measured[i] = gravity_torque[i];
          // std::cout << "t:" << gravity_torque[i] <<std::endl << "pos" << iiwa_status_.joint_position_measured[i];
        }
        lcm_.publish(kLcmGravityChannel, &grav_status);


        // PD position control
        Eigen::VectorXd position_ctrl_torque_command(kNumDof);
        Eigen::VectorXd Kp_pos_ctrl(kNumDof); // 7 joints
        // Kp_pos_ctrl << 225, 289, 144, 49, 324, 49, 49;
        // Kp_pos_ctrl << 82, 90, 94, 77, 79, 27, 37; //Mitchells gains
        Kp_pos_ctrl << 140, 295, 110, 80, 45, 40, 20;// Mitchell's gains for GPS

        Eigen::VectorXd Kd_pos_ctrl(kNumDof); // 7 joints
        // Kd_pos_ctrl << 30, 34, 24, 14, 36, 14, 14;
        // Kd_pos_ctrl << 55, 33, 22, 19, 13, 9, 5; //Mitchells gains
        Kd_pos_ctrl << 25, 33, 20, 15, 5, 2, 2; // Mitchell's gains for GPS

        // (TODOs) Add integral control (anti-windup)
        for (int joint = 0; joint < kNumJoints; joint++) {
          position_ctrl_torque_command(joint) = Kp_pos_ctrl(joint)*(joint_position_desired(joint) - iiwa_status_.joint_position_measured[joint])
                                              + Kd_pos_ctrl(joint)*(joint_velocity_desired(joint) - iiwa_status_.joint_velocity_estimated[joint]);
        }
        //Combination of ID torque control and IK position control
        torque_command += position_ctrl_torque_command;

        // -------->(For Safety) Set up iiwa position command<----------
        for (int joint = 0; joint < kNumJoints; joint++) {
          iiwa_command.joint_position[joint] = q_ref(joint);
        }

        // -------->Set up iiwa torque command<-------------
        for (int joint = 0; joint < kNumJoints; joint++) {
          iiwa_command.joint_torque[joint] = torque_command(joint);
          iiwa_command.joint_torque[joint] = std::max(-150.0, std::min(150.0, iiwa_command.joint_torque[joint]));
        }

        if(run_){
          std::cout << "I AM PUBLISHING!\n" ;
          lcm_.publish(kLcmCommandChannel, &iiwa_command);
        }
      }
    }
  }

 private:
   Eigen::VectorXd gravity_comp_no_gripper(KinematicsCache<double>& cache, const Eigen::VectorXd& vd,
       bool include_velocity_terms, const RigidBodyTree<double>& tree) const {
     cache.checkCachedKinematicsSettings(include_velocity_terms, include_velocity_terms, "gravity_comp_no_gripper");

     const bool include_acceleration_terms = true;
     int num_joints = 7;
     int kTwistSize = 6;
     unsigned int body_size_no_gripper = tree.FindBodyIndex("iiwa_link_ee") + 1; // the last arm link before gripper links, + 1 is due to additional world frame

     // Compute spatial accelerations and net wrenches that should be exerted to
     // achieve those accelerations.
     Matrix6X<double> body_accelerations(kTwistSize, body_size_no_gripper);
     Matrix6X<double> net_wrenches(kTwistSize, body_size_no_gripper);
     for (size_t i = 0; i < body_size_no_gripper; ++i) {
       const RigidBody<double>& body = *tree.bodies[i];
       if (body.has_parent_body()) {
         const RigidBody<double>& parent_body = *(body.get_parent());
         const auto& cache_element = cache.getElement(body);

         auto body_acceleration = body_accelerations.col(i);

         // Initialize body acceleration to acceleration of parent body.
         auto parent_acceleration =
             body_accelerations.col(parent_body.get_body_index());
         body_acceleration = parent_acceleration;
         // Add component due to joint acceleration.
         if (include_acceleration_terms) {
           const DrakeJoint& joint = body.getJoint();
           auto vd_joint = vd.middleRows(body.get_velocity_start_index(),
                                         joint.get_num_velocities());
           body_acceleration.noalias() +=
               cache_element.motion_subspace_in_world * vd_joint;
         }
         auto net_wrench = net_wrenches.col(i);
         const auto& body_inertia = cache_element.inertia_in_world;
         net_wrench.noalias() = body_inertia * body_acceleration;
       } else {
         drake::TwistVector<double> a_grav;
         a_grav << 0, 0, 0, 0, 0, -9.81;
         body_accelerations.col(i) = -a_grav.cast<double>();
         net_wrenches.col(i).setZero();
       }
     }

     // Do a backwards pass to compute joint wrenches from net wrenches,
     // and project the joint wrenches onto the joint's motion subspace to find the joint torque.
     auto& joint_wrenches = net_wrenches;
     const auto& joint_wrenches_const = net_wrenches;
     VectorX<double> gravity_torques(num_joints, 1);

     for (ptrdiff_t i = body_size_no_gripper - 1; i >= 0; --i) {
       RigidBody<double>& body = *tree.bodies[i];
       if (body.has_parent_body()) {
         const auto& cache_element = cache.getElement(body);
         const auto& joint = body.getJoint();
         auto joint_wrench = joint_wrenches_const.col(i);

         const auto& motion_subspace = cache_element.motion_subspace_in_world;
         auto joint_torques = gravity_torques.middleRows(body.get_velocity_start_index(), joint.get_num_velocities());
         joint_torques.noalias() = motion_subspace.transpose() * joint_wrench;

         const RigidBody<double>& parent_body = *(body.get_parent());
         auto parent_joint_wrench = joint_wrenches.col(parent_body.get_body_index());
         parent_joint_wrench += joint_wrench;
       }
     }

     return gravity_torques;
   }


  void HandleStatus(const lcm::ReceiveBuffer* rbuf, const std::string& chan,
                    const lcmt_iiwa_status* status) {
    iiwa_status_ = *status;
  }
  void HandleCancelPlan(const lcm::ReceiveBuffer* rbuf, const std::string& chan,
                    const lcmt_iiwa_status* status) {
    std::cout << "Plan Cancel Command Recieved!" << std::endl;
    run_ = false;
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
    run_ = true;
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
      drake::GetDrakePath() + "/examples/kuka_iiwa_arm/urdf/iiwa14_fixed_gripper.urdf",
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
