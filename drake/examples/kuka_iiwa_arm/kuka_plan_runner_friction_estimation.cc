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
#include <cmath>

#include <stdio.h>
#include <fstream>
#include <string>
#include <algorithm>
#include <list>
#include <iostream>

#include "robotlocomotion/robot_plan_t.hpp"

#include "drake/common/drake_assert.h"
#include "drake/common/drake_path.h"
#include "drake/common/trajectories/piecewise_polynomial.h"
#include "drake/common/trajectories/piecewise_polynomial_trajectory.h"
#include "drake/examples/kuka_iiwa_arm/iiwa_common.h"
#include "drake/multibody/parsers/urdf_parser.h"
#include "drake/multibody/rigid_body_tree.h"

#include "drake/lcmt_iiwa_status.hpp"
#include "drake/lcmt_robot_controller_reference.hpp"
#define KUKA_DATA_DIR "/home/yezhao/kuka-dev-estimation/drake/drake/examples/kuka_iiwa_arm/experiment_data/friction_model/"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using drake::Vector1d;
using Eigen::Vector2d;
using Eigen::Vector3d;

#define PI 3.1415926
static std::list< const char*> gs_fileName;
static std::list< std::string > gs_fileName_string;

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

    Eigen::VectorXd qd_meas_previous(kNumJoints); // 7DOF joint velocity at previous time sample
    qd_meas_previous.setZero();

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

        // double joint_vel = 0.6;
        // double joint_max_position = PI/6; 
        // if (joint_pos >= joint_max_position){
        //   traj_time_init_s = cur_traj_time_s;
        //   vel_sign = -1;
        //   std::cout << "hit joint pos limit" << std::endl;
        //   joint_pos_init = joint_pos - joint_offset;
        // }else if (joint_pos <= -joint_max_position){
        //   traj_time_init_s = cur_traj_time_s;
        //   vel_sign = 1;
        //   std::cout << "hit joint neg limit" << std::endl;
        //   joint_pos_init = joint_pos + joint_offset;
        // }
        // joint_pos = joint_pos_init + (double)vel_sign*joint_vel*(cur_traj_time_s - traj_time_init_s);
        
        // Eigen::VectorXd q_ref(kNumJoints);
        // Eigen::VectorXd qd_ref(kNumJoints);
        // Eigen::VectorXd qdd_ref(kNumJoints);
        // q_ref << 0,0,joint_pos,0,0,0,0;
        // qd_ref << 0,0,joint_vel,0,0,0,0;
        // qdd_ref << 0,0,0,0,0,0,0;

        Eigen::VectorXd q_meas(kNumJoints);
        Eigen::VectorXd qd_meas(kNumJoints);
        Eigen::VectorXd qdd_meas(kNumJoints);
        Eigen::VectorXd torque_meas(kNumJoints);
        for(int joint = 0; joint < kNumJoints; joint++){
          q_meas(joint) = iiwa_status_.joint_position_measured[joint];
          qd_meas(joint) = iiwa_status_.joint_velocity_estimated[joint];
          torque_meas(joint) = iiwa_status_.joint_torque_measured[joint];
        }

        //derive joint accelerations
        for (int joint = 0; joint < kNumJoints; joint++){
            qdd_meas(joint) = (qd_meas(joint) - qd_meas_previous(joint))/0.005;
        }
        qd_meas_previous = qd_meas;

        saveVector(q_meas, "joint_position_measured");
        saveVector(qd_meas, "joint_velocity_measured");
        saveVector(qdd_meas, "joint_acceleration_measured");
        saveVector(torque_meas, "joint_torque_measured");
        saveValue(cur_traj_time_s, "cur_traj_time_s");

        robot_controller_reference.utime = iiwa_status_.utime;

        for(int joint = 0; joint < kNumJoints; joint++){
          robot_controller_reference.joint_position_desired[joint] = q_ref(joint);
          robot_controller_reference.joint_velocity_desired[joint] = qd_ref(joint);
          robot_controller_reference.joint_accel_desired[joint] = qdd_ref(joint);
        }

        // publish robot controller reference to kuka control runner
        lcm_.publish(kLcmControlRefChannel, &robot_controller_reference);
      }
    }
  }

  void saveVector(const Eigen::VectorXd & _vec, const char * _name){
      std::string _file_name = KUKA_DATA_DIR;
      _file_name += _name;
      _file_name += ".txt";
      clean_file(_name, _file_name);

      std::ofstream save_file;
      save_file.open(_file_name, std::fstream::app);
      for (int i(0); i < _vec.rows(); ++i){
          save_file<<_vec(i)<< "\t";
      }
      save_file<<"\n";
      save_file.flush();
      save_file.close();
  }

  void saveValue(double _value, const char * _name){
      std::string _file_name = KUKA_DATA_DIR;
      _file_name += _name;
      _file_name += ".txt";
      clean_file(_name, _file_name);

      std::ofstream save_file;
      save_file.open(_file_name, std::fstream::app);
      save_file<<_value <<"\n";
      save_file.flush();
  }

  void clean_file(const char * _file_name, std::string & _ret_file){
      std::list<std::string>::iterator iter = std::find(gs_fileName_string.begin(), gs_fileName_string.end(), _file_name);
      if(gs_fileName_string.end() == iter){
          gs_fileName_string.push_back(_file_name);
          remove(_ret_file.c_str());
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

  auto tree = std::make_unique<RigidBodyTree<double>>();
  parsers::urdf::AddModelInstanceFromUrdfFileToWorld(
      GetDrakePath() + "/examples/kuka_iiwa_arm/urdf/iiwa14_fixed_gripper.urdf",
      multibody::joints::kFixed, tree.get());

  RobotPlanRunner runner(*tree);
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