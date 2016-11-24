#include <lcm/lcm-cpp.hpp>
#include <memory>

#include "drake/lcmt_iiwa_command.hpp"
#include "drake/lcmt_iiwa_status.hpp"

#include "drake/examples/kuka_iiwa_arm/iiwa_status.h"
#include "drake/common/drake_path.h"

namespace drake{
namespace examples{
namespace kuka_iiwa_arm{

/** Handles sending and receiving of LCM messages
 *
 */
class KukaMessageHandler {
  public:
    static const char* kLcmDarkCommandChannel;
    static const int kNumJoints;
    /** Constructor
     */
    KukaMessageHandler(){
      lcm_ = std::make_shared<lcm::LCM>();
      // subscribe the handler to IIWA_STATUS
      lcm_->subscribe(IiwaStatus<double>::channel(),
                      &KukaMessageHandler::handleMessage,
                      this);
      printf("Listening to channel: %s\n",IiwaStatus<double>::channel().c_str());
      // record the time of initialization
      time(&last_received_time_);
    }

    /** Checks if it has been more than 3 seconds since the last message
     */
    bool hasTimedOut(){
      // check if it has been more than 0.5 seconds since the last message
      return difftime(time(NULL),last_received_time_) > 3;
    }

    /** Checks for new commands from LCM
     *  This must be called in the event loop in order for new messages to come in
     */
    bool handle(){
      // wait 10ms before moving releasing control of the process
      return lcm_->handleTimeout(10);
    }

    /** Waits for the first message to come in up to a timeout
     *  Returns true if a message came in successfully and false if it failed
     */
    bool waitForFirstMessage(){
      bool success = true;
      while(!this->received_first_message_ && success){
        this->handle();
        success = !(this->hasTimedOut());
      }
      return success;
    }

    /** Prints the positions of the current state, useful for debugging
     */
    void printPositions(){
      int numJoints = iiwa_status_.num_joints;
      for (int i=0; i<numJoints; i++){
        printf("Joint %i measured at: %f\n", i, iiwa_status_.joint_position_measured[i]); 
      }
    }

    /** Returns a vector of the current positions
     */
    Eigen::VectorXd getPosition(){
      Eigen::VectorXd state(kNumJoints);
      for (int i=0; i<kNumJoints; i++){
        state[i] = iiwa_status_.joint_position_measured[i];
      }
      return state;
    }

    /** Publishes a raw lcmt_iiwa_status_command
     */
    void publish(lcmt_iiwa_command iiwa_command){
      lcm_->publish(kLcmDarkCommandChannel, &iiwa_command);
    }

    /** Publishes positions and torques from Eigen vector objects
     * Converts Eigen::VectorXd to lcm objects and publishes them
     */
    void publish(Eigen::VectorXd positions, Eigen::VectorXd torques){
      // TODO: check that the vectors are the right length and throw an exeption if they aren't

      lcmt_iiwa_command lcm_command = createCommand();

      for (int i=0; i<kNumJoints; i++){
        lcm_command.joint_position[i] = positions[i];
        lcm_command.joint_torque[i] = torques[i];
      }

      publish(lcm_command);
    }

    /** Publishes torques from an Eigen vector object
     * Converts Eigen::VectorXd to an lcm object and publishes it_interval
     * Sets positions to zero
     */
    void publishTorques(Eigen::VectorXd torques){
      // TODO: check that the vectors are the right length and throw an exeption if they aren't
      publish(Eigen::VectorXd::Zero(kNumJoints), torques);
    }

  private:
    lcmt_iiwa_status iiwa_status_;
    time_t last_received_time_;
    bool received_first_message_ = false;
    std::shared_ptr<lcm::LCM> lcm_;

    void handleMessage(const lcm::ReceiveBuffer* rbuf,
                       const std::string& chan,
                       const lcmt_iiwa_status* status){
      iiwa_status_ = *status;
      time(&last_received_time_);
      received_first_message_ = true;
    }

    lcmt_iiwa_command createCommand(){
      lcmt_iiwa_command lcm_command;
      // echo back the same timestamp
      lcm_command.timestamp = iiwa_status_.timestamp;
      lcm_command.num_joints = kNumJoints;
      lcm_command.num_torques = kNumJoints;
      lcm_command.joint_position.resize(kNumJoints, 0.);
      lcm_command.joint_torque.resize(kNumJoints, 0.);

      return lcm_command;
    }
};

const char* KukaMessageHandler::kLcmDarkCommandChannel = "IIWA_DARK_COMMAND";
const int KukaMessageHandler::kNumJoints = 7;

} // kuka_iiwa_arm
} // examples
} // drake