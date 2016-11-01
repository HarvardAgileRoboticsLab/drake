#include <gflags/gflags.h>

#include <lcm/lcm-cpp.hpp>

#include "drake/lcmt_iiwa_command.hpp"
#include "drake/lcmt_iiwa_status.hpp"

#include "drake/common/drake_assert.h"
#include "drake/common/drake_path.h"
#include "drake/common/polynomial.h"
#include "drake/common/text_logging_gflags.h"
//#include "drake/examples/kuka_iiwa_arm/iiwa_simulation.h"
#include "drake/examples/kuka_iiwa_arm/iiwa_status.h"
#include "drake/systems/LCMSystem.h"
#include "drake/systems/LinearSystem.h"
#include "drake/systems/Simulation.h"
#include "drake/systems/cascade_system.h"
#include "drake/systems/gravity_compensated_system.h"
#include "drake/systems/plants/BotVisualizer.h"

#include "drake/systems/plants/IKoptions.h"
#include "drake/systems/plants/RigidBodyIK.h"
#include "drake/systems/plants/RigidBodyTree.h"
#include "drake/systems/plants/constraint/RigidBodyConstraint.h"
#include "drake/systems/plants/joints/floating_base_types.h"
#include "drake/systems/trajectories/PiecewisePolynomial.h"
#include "drake/systems/vector.h"

using drake::AffineSystem;
using drake::BotVisualizer;
using drake::GravityCompensatedSystem;
using drake::RigidBodySystem;
using drake::SimulationOptions;
using Eigen::VectorXd;
using Eigen::MatrixXd;

namespace drake {
namespace examples {
namespace kuka_iiwa_arm {
namespace {

const char* kLcmCommandChannel = "IIWA_COMMAND";

DEFINE_double(duration, 0.75, "Simulation duration");
DEFINE_double(magnitude, 1.75, "Joint 5 Input torque magnitude");

class TrajectoryRunner {
 public:
  TrajectoryRunner(std::shared_ptr<lcm::LCM> lcm, int nT, const double* t,
                   const Eigen::MatrixXd& traj)
      : lcm_(lcm), nT_(nT), t_(t), traj_(traj) {
    lcm_->subscribe(IiwaStatus<double>::channel(),
                    &TrajectoryRunner::HandleStatus, this);
    DRAKE_ASSERT(traj_.cols() == nT);
  }

  void Run() {

    // -------->To be updated for a new motion planner<-------------
    typedef PiecewisePolynomial<double> PPType;
    typedef PPType::PolynomialType PPPoly;
    typedef PPType::PolynomialMatrix PPMatrix;
    std::vector<PPMatrix> polys;
    std::vector<double> times;

    // For each timestep, create a PolynomialMatrix for each joint
    // position.  Each column of traj_ represents a particular time,
    // and the rows of that column contain values for each joint
    // coordinate.
    for (int i = 0; i < nT_; i++) {
      PPMatrix poly_matrix(traj_.rows(), 1);
      const auto traj_now = traj_.col(i);

      // Produce interpolating polynomials for each joint coordinate.
      if (i != nT_ - 1) {
        for (int row = 0; row < traj_.rows(); row++) {
          Eigen::Vector2d coeffs(0, 0);
          coeffs[0] = traj_now(row);

          // Sets the coefficients for a linear polynomial within the interval
          // of time t_[i] < t < t_[i+1] so that the entire piecewise polynomial
          // is continuous at the time instances t_[i].
          // PiecewisePolynomial<T>::value(T t) clamps t to be between t_[0] and
          // t_[nT-1] so that for t > t_[nT-1] the piecewise polynomial always
          // evaluates to the last trajectory instance, traj_.col(nT_-1).
          coeffs[1] = (traj_(row, i + 1) - coeffs[0]) / (t_[i + 1] - t_[i]);
          poly_matrix(row) = PPPoly(coeffs);
        }
        polys.push_back(poly_matrix);
      }
      times.push_back(t_[i]);
    }

    PPType pp_traj(polys, times);

    bool time_initialized = false;
    int64_t start_time_ms = -1;
    int64_t cur_time_ms = -1;
    const int64_t end_time_offset_ms = (t_[nT_ - 1] * 1e3);
    DRAKE_ASSERT(end_time_offset_ms > 0);

    lcmt_iiwa_command iiwa_command;
    iiwa_command.num_joints = kNumJoints;
    iiwa_command.joint_position.resize(kNumJoints, 0.);
    iiwa_command.num_torques = 0;
    iiwa_command.joint_torque.resize(kNumJoints, 0.);

    while (!time_initialized ||
           cur_time_ms < (start_time_ms + end_time_offset_ms)) {
      // The argument to handleTimeout is in msec, and should be
      // safely bigger than e.g. a 200Hz input rate.
      int handled  = lcm_->handleTimeout(10);
      if (handled <= 0) {
        std::cerr << "Failed to receive LCM status." << std::endl;
        return;
      }

      if (!time_initialized) {
        start_time_ms = iiwa_status_.timestamp;
        time_initialized = true;
      }
      cur_time_ms = iiwa_status_.timestamp;

      const double cur_traj_time_s =
          static_cast<double>(cur_time_ms - start_time_ms) / 1e3;
      const auto desired_next = pp_traj.value(cur_traj_time_s);

      iiwa_command.timestamp = iiwa_status_.timestamp;

      // -------->Set up iiwa command<-------------
      // This is totally arbitrary.  There's no good reason to
      // implement this as a maximum delta to submit per tick.  What
      // we actually need is something like a proper
      // planner/interpolater which spreads the motion out over the
      // entire duration from current_t to next_t, and commands the
      // next position taking into account the velocity of the joints
      // and the distance remaining.
      const double max_joint_delta = 0.1;
      for (int joint = 0; joint < kNumJoints; joint++) {
        double joint_delta =
            desired_next(joint) - iiwa_status_.joint_position_measured[joint];
        joint_delta = std::max(-max_joint_delta,
                               std::min(max_joint_delta, joint_delta));
        iiwa_command.joint_position[joint] =
            iiwa_status_.joint_position_measured[joint] + joint_delta;
      }

      lcm_->publish(kLcmCommandChannel, &iiwa_command);
    }
  }

 private:
  void HandleStatus(const lcm::ReceiveBuffer* rbuf, const std::string& chan,
                    const lcmt_iiwa_status* status) {
    iiwa_status_ = *status;
  }

  static const int kNumJoints = 7;
  std::shared_ptr<lcm::LCM> lcm_;
  const int nT_;
  const double* t_;
  const Eigen::MatrixXd& traj_;
  lcmt_iiwa_status iiwa_status_;
};


// TODO(naveenoid) : Combine common code with
// run_kuka_iiwa_gravity_compensated_position_control into a class
// with a common method.

int main(int argc, char* argv[]) {
  std::shared_ptr<RigidBodySystem> iiwa_system = CreateKukaIiwaSystem();

  double kDuration = 0.75;
  double kInputTorqueMagnitude = 1.75;

  gflags::ParseCommandLineFlags(&argc, &argv, true);
  logging::HandleSpdlogGflags();
  kDuration = FLAGS_duration;
  kInputTorqueMagnitude = FLAGS_magnitude;

  const int kNumDof = 7;

  // Applies a small input torque on 5th joint.
  VectorXd input_torque_vector = VectorXd::Zero(kNumDof);
  input_torque_vector(4) = kInputTorqueMagnitude;

  // The input torque is generated from a constant output AffineSystem.
  // The individual matrices of the AffineSystem are all set to zero barring
  // the initial output y0 which is then set to the dimension of the inputs
  // to the IIWA System. For more details please see :
  // http://drake.mit.edu/doxygen_cxx/classdrake_1_1_affine_system.html
  auto input_torque = std::make_shared<
      AffineSystem<NullVector, NullVector, RigidBodySystem::StateVector>>(
      MatrixXd::Zero(0, 0), MatrixXd::Zero(0, 0), VectorXd::Zero(0),
      MatrixXd::Zero(kNumDof, 0), MatrixXd::Zero(kNumDof, 0),
      input_torque_vector);

  auto lcm = std::make_shared<lcm::LCM>();
  auto visualizer = CreateKukaIiwaVisualizer(iiwa_system, lcm);

  auto controlled_robot =
      std::allocate_shared<GravityCompensatedSystem<RigidBodySystem>>(
          Eigen::aligned_allocator<GravityCompensatedSystem<RigidBodySystem>>(),
          iiwa_system);

  //auto sys = cascade(cascade(input_torque, controlled_robot), visualizer);

  // Obtains an initial state of the simulation.
  //VectorXd x0 = GenerateArbitraryIiwaInitialState();
  // Specifies the start time of the simulation.
  //const double kStartTime = 0;

  //SimulationOptions options = SetupSimulation();
  //simulate(*sys.get(), kStartTime, kDuration, x0, options);

  // -------->Inverse Dynamics To be Implement<-------------
  inverseKinPointwise(&tree, kNumTimesteps, t, q0, q0, constraint_array.size(),
                      constraint_array.data(), ikoptions, &q_sol, info,
                      &infeasible_constraint);

  bool info_good = true;
  for (int i = 0; i < kNumTimesteps; ++i) {
    printf("INFO[%d] = %d ", i, info[i]);
    if (info[i] != 1) {
      info_good = false;
    }
  }
  printf("\n");

  if (!info_good) {
    std::cerr << "Solution failed, not sending." << std::endl;
    return 1;
  }

  // Now run through the plan.
  TrajectoryRunner runner(lcm, kNumTimesteps, t, q_sol);
  runner.Run();
  return 0;
}

}  // namespace
}  // namespace kuka_iiwa_arm
}  // namespace examples
}  // namespace drake

int main(int argc, char* argv[]) {
  return drake::examples::kuka_iiwa_arm::main(argc, argv);
}
