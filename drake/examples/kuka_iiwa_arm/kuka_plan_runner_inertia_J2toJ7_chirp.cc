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
#define KUKA_DATA_DIR "/home/yezhao/kuka-dev-estimation/drake/drake/examples/kuka_iiwa_arm/experiment_data/dynamic_test/"

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

static bool forward_end1 = 0;
static bool forward_end2 = 0;
static bool forward_end3 = 0;
static bool backward_end1 = 0;
static bool backward_end2 = 0;
static bool backward_end3 = 0;
static double new_init_chirp_time_s1 = 0;
static double new_init_chirp_time_s2 = 0;
static double new_init_chirp_time_s3 = 0;
static bool set_new_init_chirp_time1 = 1;
static bool set_new_init_chirp_time2 = 1;
static bool set_new_init_chirp_time3 = 1;

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

    int num_joint_pose = 7;
    int num_of_joint = 3;
    int init_joint_index = 1;
    Eigen::VectorXd joint_pose_min(num_of_joint);
    Eigen::VectorXd joint_pose_max(num_of_joint);

    joint_pose_min << -70*PI/180, -110*PI/180, -70*PI/180;
    joint_pose_max << 70*PI/180, 110*PI/180, 70*PI/180;

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

    Eigen::VectorXd q_ref(kNumJoints);
    Eigen::VectorXd qd_ref(kNumJoints);
    Eigen::VectorXd qdd_ref(kNumJoints);
    q_ref << 0,0,0,0,0,0,0;
    qd_ref << 0,0,0,0,0,0,0;
    qdd_ref << 0,0,0,0,0,0,0;

    Eigen::VectorXd joint_pos_init(num_joint_pose);
    joint_pos_init << joint_2_pose_set(3), joint_3_pose_set(3), joint_4_pose_set(3);

    Eigen::VectorXd q_meas(kNumJoints);
    Eigen::VectorXd qd_meas(kNumJoints);
    Eigen::VectorXd qdd_meas(kNumJoints);
    Eigen::VectorXd torque_meas(kNumJoints);

    Eigen::VectorXd init_joint_vel(num_of_joint);
    init_joint_vel << -0.4, -0.4, -0.4;

    Eigen::VectorXd qd_meas_previous(kNumJoints); // 7DOF joint velocity at previous time sample
    qd_meas_previous.setZero();

    double SWEEP_FREQ_HZ_LOW = 0.001;
    Eigen::VectorXd SWEEP_FREQ_HZ_HIGH(num_of_joint);
    SWEEP_FREQ_HZ_HIGH << 1, 0.7, 2;
    double SWEEP_RATE = 0.03; //percent change of sweep range per second
    double SWEEP_SETPOINT_CURRENT_AMPLITUDE_A_DEFAULT = 0.5;
    double SWEEP_SETPOINT_MID_CURRENT_A = 0.0; 

    double Switching_amp_amps = SWEEP_SETPOINT_CURRENT_AMPLITUDE_A_DEFAULT;
    double Switching_offset_amps = SWEEP_SETPOINT_MID_CURRENT_A;
    bool End1= 0;
    bool End2= 0;
    bool End3= 0;

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

        cur_traj_time_s = static_cast<double>(cur_time_us - start_time_us) / 1e6;

        double samplePeriod_s = 1e-3;

        if(!End1){
          q_ref(init_joint_index) = getChirpSignal(cur_traj_time_s, samplePeriod_s, Switching_amp_amps, Switching_offset_amps, SWEEP_FREQ_HZ_LOW, SWEEP_FREQ_HZ_HIGH(0), SWEEP_RATE, &End1, &forward_end1, &backward_end1, &new_init_chirp_time_s1, &set_new_init_chirp_time1);
        }

        if(!End2){
          q_ref(1+init_joint_index) = getChirpSignal(cur_traj_time_s, samplePeriod_s, Switching_amp_amps, Switching_offset_amps, SWEEP_FREQ_HZ_LOW, SWEEP_FREQ_HZ_HIGH(1), SWEEP_RATE, &End2, &forward_end2, &backward_end2, &new_init_chirp_time_s2, &set_new_init_chirp_time2);
        }

        if(!End3){
          q_ref(2+init_joint_index) = getChirpSignal(cur_traj_time_s, samplePeriod_s, Switching_amp_amps, Switching_offset_amps, SWEEP_FREQ_HZ_LOW, SWEEP_FREQ_HZ_HIGH(2), SWEEP_RATE, &End3, &forward_end3, &backward_end3, &new_init_chirp_time_s3, &set_new_init_chirp_time3);
        }

        for(int i = 0;i < 3;i++)
          q_ref(i+init_joint_index) += joint_pos_init(i);
        
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

        double tolerance = 0.01;
        if ((fabs(qd_meas(1)) > tolerance) || (fabs(qd_meas(2)) > tolerance) || (fabs(qd_meas(3)) > tolerance)){
          saveVector(q_meas, "joint_position_measured");
          saveVector(qd_meas, "joint_velocity_measured");
          saveVector(qdd_meas, "joint_acceleration_measured");
          saveVector(torque_meas, "joint_torque_measured");
          saveValue(cur_traj_time_s, "cur_traj_time_s");
        }

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

    //generates a chirp signal
    double getChirpSignal(double elapsedTime, double samplePeriod, double amplitude, double offset, double lowFreq_hz, double highFreq_hz, double rate, bool *end, bool *forward_end, bool *backward_end, double *new_init_chirp_time_s, bool *set_new_init_chirp_time){
      double effective_angle, effective_switching_freq_hz, chirp_signal;

      if(elapsedTime > 0.0){

        if(*backward_end){
          *end = 1;
        }
        std::cout << "------------------" << std::endl;
        if (*forward_end){
          rate = -rate;
          if(*set_new_init_chirp_time){
            *new_init_chirp_time_s = elapsedTime;
            *set_new_init_chirp_time = 0;
          }
          effective_switching_freq_hz = highFreq_hz + rate*(elapsedTime - *new_init_chirp_time_s);
        }else{
          effective_switching_freq_hz = lowFreq_hz + rate*elapsedTime;
        }
        
        effective_angle = elapsedTime * 2 * PI * effective_switching_freq_hz;
        chirp_signal = amplitude*sin(effective_angle) + offset;
        std::cout << "effective_switching_freq_hz: " << effective_switching_freq_hz << std::endl;

        if(effective_switching_freq_hz > highFreq_hz) {
          *forward_end = 1;
        }
        if((effective_switching_freq_hz < lowFreq_hz) && *forward_end) {
          *backward_end = 1;
        }
      }
      else{
        chirp_signal = offset;
        effective_switching_freq_hz = 0.0;
      }
      return chirp_signal;
    }

  void saveVector(const Eigen::VectorXd & _vec, const char * _name){
      std::string _file_name = KUKA_DATA_DIR;
      _file_name += _name;
      _file_name += ".dat";
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
      _file_name += ".dat";
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
      GetDrakePath() + "/examples/kuka_iiwa_arm/urdf/iiwa14_estimated_params_fixed_gripper.urdf",
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