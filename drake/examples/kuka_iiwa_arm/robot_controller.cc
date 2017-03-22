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
#include "drake/gps_run_controller.hpp"
#include "drake/gps_controller_gains.hpp"


using namespace std;
using ns = chrono::nanoseconds;
using get_time = chrono::steady_clock;

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
const char* const kLcmRunControllerChannel = "GPS_RUN_CONTROLLER";
const char* const kLcmCommandChannel = "IIWA_COMMAND";

  class RobotController {
  public:

    /// tree is aliased
    explicit RobotController() {
      lcm_.subscribe(kLcmStatusChannel, &RobotController::HandleStatus, this);
      lcm_.subscribe(kLcmRunControllerChannel, &RobotController::HandleRun, this);

      X_.resize(kNumStates_);
      U_.resize(kNumJoints_);

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


        while(currentStep < kNumTimeSteps_) {

          lcmt_iiwa_command iiwa_command;
          iiwa_command.num_joints = kNumJoints_;
          iiwa_command.joint_position.resize(kNumJoints_, 0.);
          iiwa_command.num_torques = kNumJoints_;
          iiwa_command.joint_torque.resize(kNumJoints_, 0.);

          U_ = K_.row(currentStep) * (X_) + k_.row(currentStep);
          vector<double> command(U_.data(), U_.data() + U_.rows() * U_.cols());
          iiwa_command.joint_torque = command;

          chrono::steady_clock::time_point end = get_time::now();
          chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(end - start);

          if (time_span.count() > dt_) {
            std::cout << "dt: "  << time_span.count() << std::endl;
            currentStep++;
            start = get_time::now();
          }

          iiwa_command.utime = iiwa_status_.utime;
          lcm_.publish(kLcmCommandChannel, &iiwa_command);
        }

        running_ = false;
      }
    }
  private:
    void HandleStatus(const lcm::ReceiveBuffer* rbuf, const std::string& chan,
                      const lcmt_iiwa_status* status) {
      iiwa_status_ = *status;
      for(int i = 0; i < kNumJoints_;i++ ) {
        X_(i) = iiwa_status_.joint_position_measured[i];
        X_(i+kNumJoints_) = iiwa_status_.joint_velocity_estimated[i];
      }
      // std::cout  << X_ << " , " << iiwa_status_.joint_position_measured[0] <<  std::endl;
    }
    void HandleRun(const lcm::ReceiveBuffer* rbuf, const std::string& chan,
                      const gps_run_controller* cmd) {
      std::cout << "GOT A NEW PLAN" << std::endl;
      if(running_ == false) {
        K_.resize(cmd->numTimeSteps, cmd->numStates);
        k_.resize(cmd->numTimeSteps, cmd->numStates);

        for(int i = 0; i < cmd->numTimeSteps; i++) {
          for(int j = 0; i < cmd->numStates; i++) {
            K_(i,j) = cmd->K.at(i).values.at(j);
            k_(i,j) = cmd->k.at(i).values.at(j);
          }
        }

        dt_ = cmd->dt;
        kNumTimeSteps_ = cmd->numTimeSteps;
        kNumJoints_ = 7; //TODO FIX! Pass via message. 7 is for iiwa
        kNumStates_ = cmd->numStates;

        X_.resize(kNumStates_);
        U_.resize(kNumJoints_);


        running_ = true;
      } else{
        std::cout << "PLAN ALREADY RUNNING. RECIEVED REQUEST IGNORED" << std::endl;
      }
    }

    VectorXd X_;
    VectorXd U_;
    MatrixXd K_;
    MatrixXd k_;
    double dt_;
    int kNumTimeSteps_;
    lcmt_iiwa_status iiwa_status_;
    lcm::LCM lcm_;
    bool running_;
    int kNumJoints_ = 7;
    int kNumStates_ = 14;

  };



 int do_main(int argc, const char* argv[]) {
   RobotController r;
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
