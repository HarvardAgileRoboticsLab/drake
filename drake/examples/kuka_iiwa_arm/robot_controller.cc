/// @file
///
/// Description!

#include <lcm/lcm-cpp.hpp>
#include <math.h>
#include <ctime>
#include <chrono>


#include "drake/common/drake_assert.h"
#include "drake/common/drake_path.h"
#include "drake/examples/kuka_iiwa_arm/iiwa_common.h"
#include "drake/multibody/parsers/urdf_parser.h"
#include "drake/multibody/rigid_body_tree.h"
#include "drake/lcmt_iiwa_command.hpp"
#include "drake/lcmt_iiwa_status.hpp"
#include "drake/lcmt_gps_trial_command.hpp"
#include "drake/lcmt_gps_controller_gains.hpp"
#include "drake/lcmt_gps_sample_result.hpp"
#include "drake/lcmt_gps_data.hpp"
#include "drake/multibody/parsers/urdf_parser.h"


using namespace std;
using ns = chrono::nanoseconds;
using get_time = chrono::steady_clock;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

enum SampleType { //IF PROTOBUF CHANGES, THIS NEEDS TO CHANGE. NEED TO INVESTIGATE IF I CAN JUST IMPORT!
  ACTION = 0,
  JOINT_ANGLES = 1,
  JOINT_VELOCITIES = 2,
  END_EFFECTOR_POINTS = 3,
  END_EFFECTOR_POINT_VELOCITIES = 4,
  END_EFFECTOR_POINT_JACOBIANS = 5,
  END_EFFECTOR_POINT_ROT_JACOBIANS = 6,
  END_EFFECTOR_POSITIONS = 7,
  END_EFFECTOR_ROTATIONS = 8,
  END_EFFECTOR_JACOBIANS = 9,
  END_EFFECTOR_HESSIANS = 10,
  RGB_IMAGE = 11,
  DEPTH_IMAGE = 12,
  RGB_IMAGE_SIZE = 13,
  CONTEXT_IMAGE = 14,
  CONTEXT_IMAGE_SIZE = 15,
  IMAGE_FEAT = 16,
  END_EFFECTOR_POINTS_NO_TARGET = 17,
  END_EFFECTOR_POINT_VELOCITIES_NO_TARGET = 18,
  TOTAL_DATA_TYPES = 19,
  JOINT_TORQUES = 20,
};


namespace drake {
namespace examples {
namespace kuka_iiwa_arm {
namespace {


const char* const kLcmStatusChannel = "IIWA_STATUS";
const char* const kLcmRunControllerChannel = "GPS_RUN_CONTROLLER";
const char* const kLcmControllerData = "GPS_DATA_RESULT";
const char* const kLcmCommandChannel = "IIWA_COMMAND";
const char* const kURDF = "IIWA_COMMAND";
const char* const BASE_FRAME = "iiwa_link_0";
const char* const EE_FRAME = "iiwa_link_ee";



  class RobotController {
  public:

    /// tree is aliased
    explicit RobotController(const RigidBodyTree<double>& tree)
        : tree_(tree) {
      lcm_.subscribe(kLcmStatusChannel, &RobotController::HandleStatus, this);
      lcm_.subscribe(kLcmRunControllerChannel, &RobotController::HandleRun, this);

      // X_.resize(T_, kNumStates_);
      // U_.resize(kNumJoints_);

      running_ = false;

    }

    void run() {
      while(true) {
        while (0 == lcm_.handleTimeout(10)) { }

        if (running_ == false || iiwa_status_.utime == 0 ) {
          continue;
        }

        int currentStep = 0;
        chrono::steady_clock::time_point start = get_time::now();

        while(currentStep < T_) {

          lcmt_iiwa_command iiwa_command;
          iiwa_command.num_joints = kNumJoints_;
          iiwa_command.joint_position.resize(kNumJoints_, 0.);
          iiwa_command.num_torques = kNumJoints_;
          iiwa_command.joint_torque.resize(kNumJoints_, 0.);

          compute_state(currentStep);

          U_ = K_.at(currentStep) * (X_.row(currentStep)) + k_.row(currentStep);

          for (int i = 0 ; i < kNumJoints_; i++) {
            iiwa_command.joint_torque.push_back(U_(i));
          }

          chrono::steady_clock::time_point end = get_time::now();
          chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(end - start);

          if (time_span.count() > dt_) {
            std::cout << "dt: "  << time_span.count()  << std::endl;
            currentStep++;
            start = get_time::now();
            compute_obs(currentStep);
          }

          iiwa_command.utime = iiwa_status_.utime;
          lcm_.publish(kLcmCommandChannel, &iiwa_command);
        }

        if(running_ == true) {
          running_ = false;
          publish_sample();
        }
      }
    }

    void compute_obs(int t){
      int idx = 0;
      for(unsigned int i = 0; i < controller_command_.obs_datatypes.size(); i++) {
        int type = controller_command_.obs_datatypes.at(i);
        switch (type) {
          case JOINT_ANGLES: {
              obs_idx_(i,0) = idx;
              for(int j = 0; j < kNumJoints_; j++) {
                X_(t, idx) = iiwa_status_.joint_position_measured[j];
                obs_idx_(i,1) = idx;
                idx++;
              }
          } case JOINT_VELOCITIES:{
              obs_idx_(i,0) = idx;
              for(int j = 0; j < kNumJoints_; j++) {
                X_(t, idx) = iiwa_status_.joint_velocity_estimated[j];
                obs_idx_(i,1) = idx;
                idx++;
              }
          } case END_EFFECTOR_POINTS: {
              obs_idx_(i,0) = idx;
              VectorXd pos;
              VectorXd vel;
              get_points_and_vels(pos,vel);
              for(int j = 0; j < pos.size(); j++) {
                X_(t,idx) = pos(j);
                obs_idx_(i,1) = idx;
                idx++;
              }
          } case END_EFFECTOR_POINT_VELOCITIES: {
              obs_idx_(i,0) = idx;
              VectorXd pos;
              VectorXd vel;
              get_points_and_vels(pos,vel);
              for(int j = 0; j < vel.size(); j++) {
                X_(t,idx) = vel(j);
                obs_idx_(i,1) = idx;
                idx++;
              }
            }
          }
        }
    }
    void publish_sample(){
      lcmt_gps_sample_result res;
      res.trialID = controller_command_.trialID;
      res.num_data_types = controller_command_.obs_datatypes.size();


      for(unsigned int i = 0; i < controller_command_.obs_datatypes.size(); i++) {
        lcmt_gps_data data;
        data.data_type = controller_command_.obs_datatypes.at(i);
        data.T = T_;
        for(int k = 0; k < T_; k++) {
          for(int j = obs_idx_(i,0); j <= obs_idx_(i,1); j++){
            data.data.at(k).push_back(obs_(k,j));
          }
        }
        res.data_samples.push_back(data);
      }

      lcmt_gps_data data;
      data.data_type = ACTION;
      for(int k = 0; k < T_; k++) {
        for(int j =0; j < kNumJoints_; k++) {
            data.data.at(k).push_back(U_(k,j));
        }
      }
      res.data_samples.push_back(data);


      for(unsigned int i = 0; i < controller_command_.state_datatypes.size(); i++) {
        lcmt_gps_data data;
        data.data_type = controller_command_.state_datatypes.at(i);

        bool already_exists = false;
        for(unsigned int a = 0; a < res.data_samples.size(); a++) {
          if(res.data_samples.at(a).data_type == data.data_type){

          }
        }
        if (already_exists) {
          continue;
        }

        data.T = T_;
        for(int k = 0; k < T_; k++) {
          for(int j = state_idx_(i,0); j <= state_idx_(i,1); j++){
            data.data.at(k).push_back(X_(k,j));
          }
        }
        res.data_samples.push_back(data);
      }

      lcm_.publish(kLcmControllerData, &res);
    }

    void compute_state(int t) {
      int idx = 0;
      for(unsigned int i = 0; i < controller_command_.state_datatypes.size(); i++) {
        int type = controller_command_.state_datatypes.at(i);
        switch (type) {
          case JOINT_ANGLES: {
              state_idx_(i,0) = idx;
              for(int j = 0; j < kNumJoints_; j++) {
                X_(t, idx) = iiwa_status_.joint_position_measured[j];
                state_idx_(i,1) = idx;
                idx++;
              }
          } case JOINT_VELOCITIES:{
            state_idx_(i,0) = idx;
              for(int j = 0; j < kNumJoints_; j++) {
                X_(t, idx) = iiwa_status_.joint_velocity_estimated[j];
                state_idx_(i,1) = idx;
                idx++;
              }
          } case END_EFFECTOR_POINTS: {
            state_idx_(i,0) = idx;
              VectorXd pos;
              VectorXd vel;
              get_points_and_vels(pos,vel);
              for(int j = 0; j < pos.size(); j++) {
                X_(t,idx) = pos(j);
                state_idx_(i,1) = idx;
                idx++;
              }
          } case END_EFFECTOR_POINT_VELOCITIES: {
            state_idx_(i,0) = idx;
              VectorXd pos;
              VectorXd vel;
              get_points_and_vels(pos,vel);
              for(int j = 0; j < vel.size(); j++) {
                X_(t,idx) = vel(j);
                state_idx_(i,1) = idx;
                idx++;
              }
            }
          }
        }
      // sample_
    }

    MatrixXd get_ee_points_in_base_frame(VectorXd ee_pos,  MatrixXd ee_rot){
      MatrixXd ee_points(kNumEEPoints_,3);

      for(int i = 0; i < ee_point_offsets_.rows(); i++) {
        ee_points.row(i) = ee_rot * ee_point_offsets_.row(i).transpose() + ee_pos.transpose();
      }

      return ee_points;
    }

    void get_ee_pos_and_rot(double * q_in, VectorXd translation, VectorXd quat, MatrixXd rotation){
      Eigen::Map<Eigen::VectorXd> q(&q_in[0], kNumJoints_); //TODO IS q the right size???? I don't think it is.
      KinematicsCache<double> cache = tree_.doKinematics(q);

      auto transform = tree_.relativeTransform(
          cache,
          tree_.FindBodyIndex(BASE_FRAME),
          tree_.FindBodyIndex(EE_FRAME));
      rotation = transform.rotation();
      translation = transform.translation();
      quat = drake::math::rotmat2quat(transform.linear());
    }

    void get_points_and_vels(VectorXd postions, VectorXd velocities){
      VectorXd posNow(3);
      VectorXd posLast(3);
      VectorXd quatNow(4);
      VectorXd quatLast(4);
      MatrixXd rotNow;
      MatrixXd rotLast;

      double qNow[kNumJoints_];
      double qLast[kNumJoints_];

      for(int i = 0; i < kNumJoints_;i++) {
        qNow[i] = iiwa_status_.joint_position_measured.at(i);
        qLast[i] = iiwa_status_old_.joint_position_measured.at(i); //Todo make this is set somewhere
      }

      get_ee_pos_and_rot(qNow, posNow, quatNow, rotNow);
      get_ee_pos_and_rot(qLast, posLast, quatLast, rotLast);

      MatrixXd pointsNow = get_ee_points_in_base_frame(posNow, rotNow);
      MatrixXd pointsLast = get_ee_points_in_base_frame(posLast, rotLast);

      MatrixXd vTmp = (pointsNow - pointsLast) / dt_;
      velocities = Eigen::Map<VectorXd>(vTmp.data(), vTmp.cols()*vTmp.rows());

      MatrixXd pTmp = pointsNow - ee_targets_;
      postions = Eigen::Map<VectorXd>(pTmp.data(), pTmp.cols()*pTmp.rows());
    }


  private:
    void HandleStatus(const lcm::ReceiveBuffer* rbuf, const std::string& chan,
                      const lcmt_iiwa_status* status) {
      iiwa_status_ = *status;
      for(int i = 0; i < kNumJoints_;i++ ) {
        // X_(i) = iiwa_status_.joint_position_measured[i];
        // X_(i+kNumJoints_) = iiwa_status_.joint_velocity_estimated[i];
      }
      // std::cout  << X_ << " , " << iiwa_status_.joint_position_measured[0] <<  std::endl;
    }
    void HandleRun(const lcm::ReceiveBuffer* rbuf, const std::string& chan,
                      const lcmt_gps_trial_command* cmd) {

      std::cout << "GOT A NEW PLAN" << std::endl;
      if(running_ == false) {

        controller_command_ = *cmd;
        T_ = cmd->T;
        dt_ = cmd->dt;
        dX_ = cmd->dX;
        dU_ = cmd->dU;
        dObs_ = cmd->dObs;
        kNumEEPoints_ = cmd->num_ee_points;
        ee_point_offsets_.resize(kNumEEPoints_ ,3);

        for(int i = 0; i < kNumEEPoints_ / 3.0; i++) {
          // std::cout << i << ", " << cmd->ee_points.size()  << ", " <<  kNumEEPoints_ / 3.0<<  std::endl;
          ee_point_offsets_(i,0) = cmd->ee_points.at((i*3));
          ee_point_offsets_(i,1) = cmd->ee_points.at((i*3)+1);
          ee_point_offsets_(i,2) = cmd->ee_points.at((i*3)+2);
        }

        ee_targets_.resize(kNumEEPoints_,3);
        for(int i = 0; i < kNumEEPoints_ / 3.0; i++) {
          ee_point_offsets_(i,0) = cmd->ee_points_tgt.at((i*3));
          ee_point_offsets_(i,1) = cmd->ee_points_tgt.at((i*3)+1);
          ee_point_offsets_(i,2) = cmd->ee_points_tgt.at((i*3)+2);
        }

        for(int i = 0; i < T_; i++) {
          K_.push_back(Eigen::MatrixXd(dU_,dX_));
        }

        k_.resize(T_, dX_);
        X_.resize(T_,dX_);
        obs_.resize(T_, dObs_);
        U_.resize(T_,dU_);
        obs_idx_.resize(cmd->obs_datatypes.size(),2);
        state_idx_.resize(cmd->state_datatypes.size(),2);


        for(int i = 0; i < T_; i++) {
          for(int j = 0; j < dU_; j++) {
            k_(i,j) = cmd->k.at(i).k.at(j);
          }
        }

        for(int i = 0; i < T_; i++) {
          for(int j = 0; j < dU_; j++) {
            for(int k = 0; k < dX_; k++) {
              K_.at(i)(j,k) = cmd->K.at(i).K.at(j).at(k);
            }
          }
        }


        kNumJoints_ = 7; //TODO FIX! Pass via message. 7 is for iiwa
        // X_.resize(dX_);
        // U_.resize(kNumJoints_);
        running_ = true;

      } else{
        std::cout << "PLAN ALREADY RUNNING. RECIEVED REQUEST IGNORED" << std::endl;
      }
    }

    std::vector<std::vector<double>> state_;

    MatrixXd X_;
    MatrixXd obs_;
    MatrixXd obs_idx_;
    MatrixXd state_idx_;

    VectorXd U_;
    std::vector<MatrixXd> K_;
    MatrixXd k_;
    MatrixXd ee_point_offsets_;
    MatrixXd ee_targets_;
    double dt_;
    int dX_;
    int dU_;
    int dObs_;

    lcmt_iiwa_status iiwa_status_;
    lcmt_iiwa_status iiwa_status_old_;
    lcm::LCM lcm_;
    bool running_;
    int kNumJoints_ = 7;
    int kNumStates_ = 14;
    int kNumEEPoints_;
    int T_;
    int trialID_ = -1;
    lcmt_gps_trial_command controller_command_;
    lcmt_gps_sample_result sample_;
    const RigidBodyTree<double>& tree_;

  };



 int do_main(int argc, const char* argv[]) {

   auto tree = std::make_unique<RigidBodyTree<double>>();
   parsers::urdf::AddModelInstanceFromUrdfFileToWorld(
       GetDrakePath() + "/examples/kuka_iiwa_arm/urdf/iiwa14_fixed_gripper.urdf",
       multibody::joints::kFixed, tree.get());


   RobotController r(*tree);
   r.run();

   return 0;
 }

 }  // namespace
 }  // namespace kuka_iiwa_arm
 }  // namespace examples
 }  // namespace drake

 int main(int argc, const char* argv[]) {
   return drake::examples::kuka_iiwa_arm::do_main(argc, argv);
 }
