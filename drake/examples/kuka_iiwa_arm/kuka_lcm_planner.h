#include <lcm/lcm-cpp.hpp>
#include <memory>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <time.h>


#include "robotlocomotion/robot_plan_t.hpp"
#include "bot_core/twist_t.hpp"
#include "bot_core/position_3d_t.hpp"
#include "bot_core/quaternion_t.hpp"
#include "bot_core/vector_3d_t.hpp"
#include "bot_core/force_torque_t.hpp"
#include "drake/lcmt_generic_planner_request.hpp"
#include "drake/lcmt_matlab_plan_request.hpp"
#include "drake/lcmt_matlab_plan_response.hpp"


#include "drake/multibody/ik_options.h"
#include "drake/multibody/rigid_body_ik.h"
#include "drake/multibody/rigid_body_tree.h"
#include "drake/common/drake_path.h"

#include "rigid_body_plan_parser.h"
// #include "kuka_dircol_params.h"

namespace drake{
namespace examples{
namespace kuka_iiwa_arm{
namespace{
bot_core::robot_state_t lcmRobotState(double, Eigen::VectorXd, RigidBodyTree<double>*);

class KukaPlanner{
  public:
    static const char* IK_REQUEST_CHANNEL;
    static const char* PLAN_REQUEST_CHANNEL;
    static const char* IK_RESPONSE_CHANNEL;
    static const char* PLAN_RESPONSE_CHANNEL;

    KukaPlanner(RigidBodyTree<double>* kuka,
                std::shared_ptr<lcm::LCM> lcm){
      lcm_ = lcm;
      kuka_ = kuka;

      lcm_->subscribe(IK_REQUEST_CHANNEL,
                      &KukaPlanner::handleRequest,
                      this);
      lcm_->subscribe(PLAN_REQUEST_CHANNEL,
                      &KukaPlanner::handleRequest,
                      this);
    }

    void run(){
      std::cout << "Starting the main loop\n";
      while(true){
        lcm_->handleTimeout(1000);
      }
    }

    void handleRequest(const lcm::ReceiveBuffer* rbuf,
                       const std::string& chan,
                       const lcmt_generic_planner_request* status){
      std::cout << "Messge from " << chan << std::endl;
      std::cout << "Handling request in KukaPlanner, passing control to appropriate function\n";
      
      if (chan == IK_REQUEST_CHANNEL){
        this->handleIkRequest(status);
      }else if (chan == PLAN_REQUEST_CHANNEL){
        this->handlePlanRequest(status);
      }
    }

  protected:
    RigidBodyTree<double>* kuka_;
    std::shared_ptr<lcm::LCM> lcm_;

    virtual void handleIkRequest(const lcmt_generic_planner_request*) = 0;

    virtual void handlePlanRequest(const lcmt_generic_planner_request*) = 0;

    void publishIkResponse(Eigen::VectorXd* pose, int info){
      // build the state
      robotlocomotion::robot_plan_t plan;
      plan.utime = 0;
      plan.robot_name = "iiwa"; 
      plan.num_states = 1;
      std::cout << "tag1" << std::endl;
      Eigen::VectorXd state(kuka_->get_num_positions() + kuka_->get_num_velocities());
      state.head(kuka_->get_num_positions()) = *pose;
      state.tail(kuka_->get_num_velocities()) = Eigen::VectorXd::Zero(kuka_->get_num_velocities());
      std::cout << "tag2" << std::endl;
      
      plan.plan.push_back(lcmRobotState(0.0,state,kuka_));
      std::cout << "tag3" << std::endl;

      plan.plan_info.push_back(info);
      // publish
      std::cout << "tag4" << std::endl;
      lcm_->publish(IK_RESPONSE_CHANNEL, &plan);
      std::cout << "tag5" << std::endl;
    }

    void publishPlanResponse(Eigen::VectorXd* t_vec, Eigen::MatrixXd* trajectory, std::vector<int> info){
      robotlocomotion::robot_plan_t plan;
      plan.utime=0;
      plan.robot_name="iiwa";
      plan.num_states = t_vec->size();
      
      for (int i=0; i<plan.num_states; i++){
        plan.plan.push_back(lcmRobotState((*t_vec)[i], trajectory->col(i), kuka_));
        plan.plan_info.push_back(info[i]);
      }
      plan.num_grasp_transitions = 0;

      plan.left_arm_control_type = plan.NONE;
      plan.left_leg_control_type = plan.NONE;
      plan.right_arm_control_type = plan.NONE;
      plan.right_leg_control_type = plan.NONE;

      plan.num_bytes = 0;

      // TODO: build the message from the time and trajectory
      std::cout << "about to publish plan" << std::endl;
      lcm_->publish(PLAN_RESPONSE_CHANNEL, &plan);
    }

};

const char* KukaPlanner::IK_REQUEST_CHANNEL = "IK_REQUEST";
const char* KukaPlanner::PLAN_REQUEST_CHANNEL = "PLANNER_REQUEST";
const char* KukaPlanner::IK_RESPONSE_CHANNEL = "CANDIDATE_MANIP_IKPLAN";
const char* KukaPlanner::PLAN_RESPONSE_CHANNEL = "CANDIDATE_MANIP_PLAN";

class KukaIkPlanner : public KukaPlanner{
  public:
    KukaIkPlanner(RigidBodyTree<double>* kuka,
                 std::shared_ptr<lcm::LCM> lcm) : KukaPlanner(kuka, lcm){}

  protected:
    virtual void handleIkRequest(const lcmt_generic_planner_request* status) {
      std::cout << "Handling ik request from KukaIkPlanner\n";

      // unpack the LCM message
      auto constraints = parse_constraints(status->constraints, kuka_);
      std::vector<RigidBodyConstraint*> constraint_ptrs(constraints.size());
      for (unsigned int i=0; i<constraints.size(); i++){
        constraint_ptrs[i] = constraints[i].get();
      }
      auto poses = parse_json_object(status->poses);
      auto joint_names = parse_json_list(status->joint_names);
      auto seed_pose_full = parse_json_list(poses[status->seed_pose]);
      auto nominal_pose_full = parse_json_list(poses[status->nominal_pose]);
      Eigen::VectorXd seed_pose(kuka_->get_num_positions());
      Eigen::VectorXd nominal_pose(kuka_->get_num_positions());
      auto pos_idx_map = kuka_->computePositionNameToIndexMap();
      for (unsigned int i=0; i<joint_names.size(); i++){
        auto joint = joint_names[i]; 
        // ignore all of the positions associated with the floating base
        if (contains(joint,"base"))
          continue;
        int idx = pos_idx_map[joint];
        seed_pose[idx] = parse_json_double(seed_pose_full[i]);
        nominal_pose[idx] = parse_json_double(nominal_pose_full[i]);
      }
      // ignore the specified ikoptions for now
      // TODO: read the ikoptions from the message
      IKoptions ikoptions(kuka_);

      //// print some of the message components for debugging
      std::cout << "---------- Message Values ------------" << std::endl;
      std::cout << "Poses: " << status->poses << std::endl;
      std::cout << "Seed Pose: " << status->seed_pose << std::endl;
      std::cout << "Nominal Pose: " << status->nominal_pose << std::endl;
      std::cout << "End Pose: " << status->end_pose << std::endl;
      std::cout << "Joint Names: " << status->joint_names << std::endl;
      std::cout << "Options: " << status->options << std::endl;
      // print some of the computed values
      std::cout << "---------- Computed Values ------------" << std::endl;
      std::cout << "Seed Pose: \n" << seed_pose << std::endl;
      std::cout << "Nominal Pose: \n" << nominal_pose << std::endl;

      auto results = inverseKinSimple(kuka_, seed_pose, nominal_pose, constraint_ptrs, ikoptions);
      std::cout << "Computed IK" << std::endl;
      publishIkResponse(&(results.q_sol[0]),results.info[0]);
    }

    virtual void handlePlanRequest(const lcmt_generic_planner_request* status) {
      std::cout << "Handling plan request in KukaIKPlanner\n";  
      const int num_timesteps = 20;
      Eigen::VectorXd time_vec = Eigen::VectorXd::Zero(num_timesteps);
      Eigen::MatrixXd traj = Eigen::MatrixXd::Zero(kuka_->get_num_positions(), num_timesteps);
      std::vector<int> info(num_timesteps);
      for (unsigned int i=0; i<num_timesteps; i++){
        info[i]=1;
      }
      // TODO compute IK and publish response
      std::cout << "Publishing an empty trajectory:\n" << traj << std::endl;
      publishPlanResponse(&time_vec, &traj, info);
      std::cout << "Finished publishing\n";
    }
};


class KukaDircolPlanner : public KukaIkPlanner {
  public:
    KukaDircolPlanner(RigidBodyTree<double>* kuka,
                 std::shared_ptr<lcm::LCM> lcm) : KukaIkPlanner(kuka, lcm){}

  protected:
    // override the plan request, but not the IK request
    virtual void handlePlanRequest(const lcmt_generic_planner_request* status) {
      std::cout << "Handling plan request from KukaDircolPlanner\n";
      // unpack the LCM message: there are usually no constraints in this plan request

      auto poses = parse_json_object(status->poses);
      auto joint_names = parse_json_list(status->joint_names);
      auto start_pose_full = parse_json_list(poses[status->seed_pose]);
      auto nominal_pose_full = parse_json_list(poses[status->nominal_pose]);
      auto goal_pose_full = parse_json_list(poses[status->end_pose]);
      Eigen::VectorXd start_pose(kuka_->get_num_positions());
      Eigen::VectorXd nominal_pose(kuka_->get_num_positions());
      Eigen::VectorXd goal_pose(kuka_->get_num_positions());
      auto pos_idx_map = kuka_->computePositionNameToIndexMap();
      for (unsigned int i=0; i<joint_names.size(); i++){
        auto joint = joint_names[i]; 
        // ignore all of the positions associated with the floating base
        if (contains(joint,"base"))
          continue;
        int idx = pos_idx_map[joint];
        start_pose[idx] = parse_json_double(start_pose_full[i]);
        nominal_pose[idx] = parse_json_double(nominal_pose_full[i]);
        goal_pose[idx] = parse_json_double(goal_pose_full[i]);
      }

      
      // TODO: compute the dynamic trajectory
      
      const int num_timesteps = 20;

      // print some of the message components for debugging
      std::cout << "---------- Message Values ------------" << std::endl;
      std::cout << "Constraints" << status->constraints << std::endl;
      std::cout << "Poses: " << status->poses << std::endl;
      std::cout << "Seed Pose: " << status->seed_pose << std::endl;
      std::cout << "Nominal Pose: " << status->nominal_pose << std::endl;
      std::cout << "End Pose: " << status->end_pose << std::endl;
      std::cout << "Joint Names: " << status->joint_names << std::endl;
      std::cout << "Options: " << status->options << std::endl;
      // print some of the computed values
      std::cout << "---------- Computed Values ------------" << std::endl;
      std::cout << "Start Pose: \n" << start_pose << std::endl;
      std::cout << "Nominal Pose: \n" << nominal_pose << std::endl;
      std::cout << "Goal Pose: \n" << goal_pose << std::endl;

      // Dummy trajectory for now
      Eigen::VectorXd time_vec = Eigen::VectorXd::Zero(num_timesteps);
      Eigen::MatrixXd traj = Eigen::MatrixXd::Zero(kuka_->get_num_positions(), num_timesteps);
      std::vector<int> info(num_timesteps);
      for (int i=0; i<num_timesteps; i++){
        info[i]=1;
      }

      std::cout << "Publishing an empty trajectory:\n" << traj << std::endl;
      publishPlanResponse(&time_vec, &traj, info);
      std::cout << "Finished publishing\n";
    }
};

class KukaMatlabDircolPlanner : public KukaIkPlanner {
  public:
    static const char* MATLAB_PLAN_REQUEST_CHANNEL;
    static const char* MATLAB_PLAN_RESPONSE_CHANNEL;

    KukaMatlabDircolPlanner(RigidBodyTree<double>* kuka,
                 std::shared_ptr<lcm::LCM> lcm) : KukaIkPlanner(kuka, lcm){
      
      lcm_->subscribe(MATLAB_PLAN_RESPONSE_CHANNEL,
        &KukaMatlabDircolPlanner::handleMatlabResponse, this);
    }

    void handleMatlabResponse(const lcm::ReceiveBuffer* rbuf,
                              const std::string& chan,
                              const lcmt_matlab_plan_response* plan){

      int N = plan->num_timesteps;
      int Nx = kuka_->get_num_positions() + kuka_->get_num_velocities();
      Eigen::VectorXd t(N);
      Eigen::MatrixXd traj(Nx, N);
      std::vector<int> info;

      for (int i=0; i<N; i++){
        std::cout << "parsing time " << i << std::endl;
        t[i] = plan->time[i];
        std::cout << "parsing position " << i << std::endl;
        info.push_back(1);//plan->status;
        for (int j=0; j < Nx; j++){
          traj(j,i) = plan->state[j][i];
        }
      }

      publishPlanResponse(&t, &traj, info);
    }
  protected:
    // override the plan request, but not the IK request
    virtual void handlePlanRequest(const lcmt_generic_planner_request* status) {
      auto poses = parse_json_object(status->poses);
      auto joint_names = parse_json_list(status->joint_names);
      auto start_pose_full = parse_json_list(poses[status->seed_pose]);
      auto nominal_pose_full = parse_json_list(poses[status->nominal_pose]);
      auto goal_pose_full = parse_json_list(poses[status->end_pose]);
      Eigen::VectorXd start_pose(kuka_->get_num_positions());
      Eigen::VectorXd nominal_pose(kuka_->get_num_positions());
      Eigen::VectorXd goal_pose(kuka_->get_num_positions());
      auto pos_idx_map = kuka_->computePositionNameToIndexMap();
      for (unsigned int i=0; i<joint_names.size(); i++){
        auto joint = joint_names[i]; 
        // ignore all of the positions associated with the floating base
        if (contains(joint,"base"))
          continue;
        int idx = pos_idx_map[joint];
        start_pose[idx] = parse_json_double(start_pose_full[i]);
        nominal_pose[idx] = parse_json_double(nominal_pose_full[i]);
        goal_pose[idx] = parse_json_double(goal_pose_full[i]);
      }

      std::cout << "Start Pose" << poses[status->seed_pose] << std::endl;
      std::cout << "End Pose" << poses[status->end_pose] << std::endl;
      lcmt_matlab_plan_request msg;
      msg.timestamp = time(NULL);
      msg.state_size = kuka_->get_num_positions()*2;
      // set position
      for (int i=0; i<kuka_->get_num_positions(); i++){
        msg.start_state.push_back(start_pose[i]);
        msg.goal_state.push_back(goal_pose[i]);
      }
      // set velocity
      for (int i=0; i<kuka_->get_num_positions(); i++){
        msg.start_state.push_back(0);
        msg.goal_state.push_back(0);
      }

      // send request
      lcm_->publish(MATLAB_PLAN_REQUEST_CHANNEL, &msg);

    }
};

const char* KukaMatlabDircolPlanner::MATLAB_PLAN_REQUEST_CHANNEL = "MATLAB_KUKA_DIRCOL_REQUEST";
const char* KukaMatlabDircolPlanner::MATLAB_PLAN_RESPONSE_CHANNEL = "MATLAB_KUKA_DIRCOL_RESPONSE";

bot_core::robot_state_t lcmRobotState(double t, Eigen::VectorXd q, RigidBodyTree<double>* tree){
  
  auto pos = q.head(tree->get_num_positions());
  auto vel = q.tail(tree->get_num_velocities());

  bot_core::robot_state_t msg;
  
  msg.utime = (int) (t*1e6);
  msg.num_joints = tree->get_num_positions();

  for (int i=0; i<msg.num_joints; i++){
    msg.joint_name.push_back(tree->get_position_name(i));
    msg.joint_position.push_back(pos[i]);
    msg.joint_velocity.push_back(vel[i]);
    msg.joint_effort.push_back(0.0);
  }

  // populate all of the unused fields to avoid seg-faults
  // pose
  bot_core::position_3d_t pose;
  bot_core::vector_3d_t translation;
  bot_core::quaternion_t rotation;

  translation.x = 0.0;
  translation.y = 0.0;
  translation.z = 0.0;

  rotation.w = 0.0;
  rotation.x = 0.0;
  rotation.y = 0.0;
  rotation.z = 0.0;

  pose.translation = translation;
  pose.rotation = rotation;
  msg.pose = pose;

  // twist
  bot_core::twist_t twist;
  bot_core::vector_3d_t linear_vel;
  bot_core::vector_3d_t angular_vel;

  linear_vel.x = 0.0;
  linear_vel.y = 0.0;
  linear_vel.z = 0.0;

  angular_vel.x = 0.0;
  angular_vel.y = 0.0;
  angular_vel.z = 0.0;

  twist.linear_velocity = linear_vel;
  twist.angular_velocity = angular_vel;

  msg.twist = twist;

  // force_torque
  bot_core::force_torque_t force_torque;
  force_torque.l_foot_force_z = 0.0;
  force_torque.l_foot_torque_x = 0.0;
  force_torque.l_foot_torque_y = 0.0;

  force_torque.r_foot_force_z = 0.0;
  force_torque.r_foot_torque_x = 0.0;
  force_torque.r_foot_torque_y = 0.0;

  for (int i=0; i<3; i++){
    force_torque.l_hand_force[i] = 0.0;  
    force_torque.l_hand_torque[i] = 0.0;
    force_torque.r_hand_force[i] = 0.0;  
    force_torque.r_hand_torque[i] = 0.0; 
  }

  msg.force_torque = force_torque;

  return msg;
}

} // anonymous
} // kuka_iiwa_arm
} // examples
} // drake