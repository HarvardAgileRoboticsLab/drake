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
#define KUKA_DATA_DIR "/home/yezhao/kuka-dev-estimation/drake/drake/examples/kuka_iiwa_arm/experiment_data/static_test/"

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

    int num_joint_pose = 4;
    int num_of_joint = 3;
    int init_joint_index = 1;
    Eigen::VectorXd joint_pose_min(num_of_joint);
    Eigen::VectorXd joint_pose_max(num_of_joint);

    joint_pose_min << -60*PI/180, -110*PI/180, -60*PI/180;
    joint_pose_max << 60*PI/180, 110*PI/180, 60*PI/180;

    Eigen::VectorXd joint_pose_inc(num_of_joint);
    joint_pose_inc << (joint_pose_max(0) - joint_pose_min(0))/(double)(num_joint_pose-1), (joint_pose_max(1) - joint_pose_min(1))/(double)(num_joint_pose-1), (joint_pose_max(2) - joint_pose_min(2))/(double)(num_joint_pose-1);

    std::cout << "joint_pose_inc" << joint_pose_inc << std::endl;

    Eigen::VectorXd joint_2_pose_set(num_joint_pose);
    Eigen::VectorXd joint_3_pose_set(num_joint_pose);
    Eigen::VectorXd joint_4_pose_set(num_joint_pose);

    for(int i = 0;i<num_joint_pose;i++){
      joint_2_pose_set(i) = joint_pose_min(0) + joint_pose_inc(0)*i;
      joint_3_pose_set(i) = joint_pose_min(1) + joint_pose_inc(1)*i;
      joint_4_pose_set(i) = joint_pose_min(2) + joint_pose_inc(2)*i;
    }

    std::cout << "joint_2_pose_set" << joint_2_pose_set << std::endl;
    std::cout << "joint_3_pose_set" << joint_3_pose_set << std::endl;
    std::cout << "joint_4_pose_set" << joint_4_pose_set << std::endl;

    int trajectory_seg_index = 1;
    int first_trajectory_seg_index_offset = 0;
    int second_trajectory_seg_index_offset = 0;
    double trajectory_seg_duration = 6;

    Eigen::VectorXd q_ref(kNumJoints);
    Eigen::VectorXd qd_ref(kNumJoints);
    Eigen::VectorXd qdd_ref(kNumJoints);
    q_ref << 0,0,0,0,0,0,0;
    qd_ref << 0,0,0,0,0,0,0;
    qdd_ref << 0,0,0,0,0,0,0;

    double traj_time_init_s = 0;
    Eigen::VectorXd joint_pos_init(num_joint_pose);
    joint_pos_init << joint_2_pose_set(0), joint_3_pose_set(0), joint_4_pose_set(0);
    bool save_initial_data_flag_ = true;
    bool ending_motion_flag_ = false;

    Eigen::VectorXd q_meas(kNumJoints);
    Eigen::VectorXd qd_meas(kNumJoints);
    Eigen::VectorXd qdd_meas(kNumJoints);
    Eigen::VectorXd torque_meas(kNumJoints);

    Eigen::VectorXd init_joint_vel(num_of_joint);
    init_joint_vel << -0.4, -0.4, -0.4;

    Eigen::VectorXd qd_meas_previous(kNumJoints); // 7DOF joint velocity at previous time sample
    qd_meas_previous.setZero();

    int first_joint_index = 0;
    int second_joint_index = 0;
    int third_joint_index = 0;

    bool first_joint_forward_index = true;
    bool second_joint_forward_index = true;
    bool third_joint_forward_index = true;

    while (true) {
      // Call lcm handle until at least one message is processed
      while (0 == lcm_.handleTimeout(10)) { }

      DRAKE_ASSERT(iiwa_status_.utime != -1);
      cur_time_us = iiwa_status_.utime;

      double cur_traj_time_s = static_cast<double>(cur_time_us - start_time_us) / 1e6;

      if (qtraj_) {
        if (plan_number_ != cur_plan_number) {
          std::cout << "Starting new plan." << std::endl;
          start_time_us = cur_time_us;
          cur_plan_number = plan_number_;
        }

        double joint_vel = 0.3;

        if (cur_traj_time_s < 10){
            for(int i = 0;i<num_of_joint;i++){
              if (abs(iiwa_status_.joint_position_measured[i+init_joint_index] - joint_pos_init(i)) > 1e-4){
                q_ref(i+init_joint_index) = init_joint_vel(i) * (cur_traj_time_s - traj_time_init_s); 
              }else{
                q_ref(i+init_joint_index) = joint_pos_init(i);
              }
            }
        }else{
          cur_traj_time_s = cur_traj_time_s - 10;
          if ((cur_traj_time_s > trajectory_seg_duration*trajectory_seg_index || save_initial_data_flag_) & !ending_motion_flag_){

            // save measured data
            for(int joint = 0; joint < kNumJoints; joint++){
              q_meas(joint) = iiwa_status_.joint_position_measured[joint];
              qd_meas(joint) = iiwa_status_.joint_velocity_estimated[joint];
              torque_meas(joint) = iiwa_status_.joint_torque_measured[joint];
            }

            //derive joint accelerations
            for (int joint = 0; joint < kNumJoints; joint++){
                qdd_meas(joint) = (qd_meas(joint) - qd_meas_previous(joint))/0.001;
            }
            qd_meas_previous = qd_meas;

            saveVector(q_meas, "joint_position_measured");
            saveVector(qd_meas, "joint_velocity_measured");
            saveVector(qdd_meas, "joint_acceleration_measured");
            saveVector(torque_meas, "joint_torque_measured");
            saveValue(cur_traj_time_s, "cur_traj_time_s");

            if (save_initial_data_flag_)
              save_initial_data_flag_ = false;
            else{
              traj_time_init_s = cur_traj_time_s;

              std::cout << "----------start----------" << std::endl;
              if (trajectory_seg_index < num_joint_pose*num_joint_pose*num_joint_pose){
                std::cout << "trajectory_seg_index:" << trajectory_seg_index << std::endl;
                first_joint_index = (int)floor((double)trajectory_seg_index/(double)(num_joint_pose*num_joint_pose));
                second_joint_index = (int)floor((double)trajectory_seg_index/(double)num_joint_pose) % num_joint_pose;
                third_joint_index = trajectory_seg_index % num_joint_pose;

                std::cout << "first_joint_forward_index: " << first_joint_forward_index << std::endl;
                std::cout << "second_joint_forward_index: " << second_joint_forward_index << std::endl;
                std::cout << "third_joint_forward_index: " << third_joint_forward_index << std::endl;

                if (trajectory_seg_index - first_trajectory_seg_index_offset > num_joint_pose*num_joint_pose){
                  if (first_joint_index == 0){
                    first_joint_forward_index = !first_joint_forward_index;
                  }
                  first_trajectory_seg_index_offset = trajectory_seg_index;
                  second_trajectory_seg_index_offset = trajectory_seg_index;
                }
                std::cout << "second_trajectory_seg_index_offset: " << second_trajectory_seg_index_offset << std::endl;
                std::cout << "second_joint_index: " << second_joint_index << std::endl;

                if (trajectory_seg_index - second_trajectory_seg_index_offset > num_joint_pose + 2){
                  if (second_joint_index == 0){
                    second_joint_forward_index = !second_joint_forward_index;
                  }
                }

                if (trajectory_seg_index >= num_joint_pose-1){    
                  if (third_joint_index == 0){
                    third_joint_forward_index = !third_joint_forward_index;
                  }
                }

                std::cout << "first_joint_forward_index: " << first_joint_forward_index << std::endl;
                std::cout << "second_joint_forward_index: " << second_joint_forward_index << std::endl;
                std::cout << "third_joint_forward_index: " << third_joint_forward_index << std::endl;

                if (first_joint_forward_index){
                  joint_pos_init(0) = joint_3_pose_set(first_joint_index);
                }else{
                  joint_pos_init(0) = joint_3_pose_set(num_joint_pose - first_joint_index - 1);
                }
                
                if (second_joint_forward_index){
                  joint_pos_init(1) = joint_3_pose_set(second_joint_index);
                }else{
                  joint_pos_init(1) = joint_3_pose_set(num_joint_pose - second_joint_index - 1);
                }

                if (third_joint_forward_index){
                  joint_pos_init(2) = joint_4_pose_set(third_joint_index);
                }else{
                  joint_pos_init(2) = joint_4_pose_set(num_joint_pose - third_joint_index - 1);
                }
                joint_pos_init(0) = joint_2_pose_set(first_joint_index);

                //joint_pos_init << joint_2_pose_set(first_joint_index), joint_3_pose_set(second_joint_index), joint_4_pose_set(third_joint_index);

                if (first_joint_forward_index)
                  std::cout << "first_joint_index: " << first_joint_index;
                else
                  std::cout << "first_joint_index: " << num_joint_pose - first_joint_index - 1;

                if (second_joint_forward_index)
                  std::cout << ", second_joint_index: " << second_joint_index;
                else
                  std::cout << ", second_joint_index: " << num_joint_pose - second_joint_index - 1;

                if (third_joint_forward_index)
                  std::cout << ", third_joint_index: " << third_joint_index << std::endl;
                else
                  std::cout << ", third_joint_index: " << num_joint_pose - third_joint_index - 1 << std::endl;
                                
                trajectory_seg_index++;
              }else{
                ending_motion_flag_ = true;
              }
            }
          }

          if (first_joint_forward_index){
            if (q_ref(init_joint_index) < joint_pos_init(0) + joint_pose_inc(0)){
              q_ref(init_joint_index) = joint_pos_init(0) + joint_vel * (cur_traj_time_s - traj_time_init_s); 
            }else{
              q_ref(init_joint_index) = joint_pos_init(0) + joint_pose_inc(0);
            } 
          }else{
            if (q_ref(init_joint_index) > joint_pos_init(0) - joint_pose_inc(0)){
              q_ref(init_joint_index) = joint_pos_init(0) - joint_vel * (cur_traj_time_s - traj_time_init_s); 
            }else{
              q_ref(1+init_joint_index) = joint_pos_init(0) - joint_pose_inc(0);
            }
          }
          
          if (second_joint_forward_index){
            if (q_ref(1+init_joint_index) < joint_pos_init(1) + joint_pose_inc(1)){
              q_ref(1+init_joint_index) = joint_pos_init(1) + joint_vel * (cur_traj_time_s - traj_time_init_s); 
            }else{
              q_ref(1+init_joint_index) = joint_pos_init(1) + joint_pose_inc(1);
            } 
          }else{
            if (q_ref(1+init_joint_index) > joint_pos_init(1) - joint_pose_inc(1)){
              q_ref(1+init_joint_index) = joint_pos_init(1) - joint_vel * (cur_traj_time_s - traj_time_init_s); 
            }else{
              q_ref(1+init_joint_index) = joint_pos_init(1) - joint_pose_inc(1);
            }
          }

          if (third_joint_forward_index){
            if (q_ref(2+init_joint_index) < joint_pos_init(2) + joint_pose_inc(2)){
              q_ref(2+init_joint_index) = joint_pos_init(2) + joint_vel * (cur_traj_time_s - traj_time_init_s); 
            }else{
              q_ref(2+init_joint_index) = joint_pos_init(2) + joint_pose_inc(2);
            }
          }else{
            if (q_ref(2+init_joint_index) > joint_pos_init(2) - joint_pose_inc(2)){
              q_ref(2+init_joint_index) = joint_pos_init(2) - joint_vel * (cur_traj_time_s - traj_time_init_s); 
            }else{
              q_ref(2+init_joint_index) = joint_pos_init(2) - joint_pose_inc(2);
            }
          }

        }

        // save measured data
        for(int joint = 0; joint < kNumJoints; joint++){
          q_meas(joint) = iiwa_status_.joint_position_measured[joint];
          qd_meas(joint) = iiwa_status_.joint_velocity_estimated[joint];
          torque_meas(joint) = iiwa_status_.joint_torque_measured[joint];
        }

        //derive joint accelerations
        for (int joint = 0; joint < kNumJoints; joint++){
            qdd_meas(joint) = (qd_meas(joint) - qd_meas_previous(joint))/0.001;
        }
        qd_meas_previous = qd_meas;

/*        saveValue(cur_traj_time_s, "cur_traj_time_s_full");
        saveVector(q_meas, "joint_position_measured_full");
        saveVector(qd_meas, "joint_velocity_measured_full");
        saveVector(torque_meas, "joint_torque_measured_full");
*/
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