/// @file
///
/// Description!

#include <lcm/lcm-cpp.hpp>
#include <math.h>

#include "drake/common/drake_assert.h"
#include "drake/common/drake_path.h"
#include "drake/examples/kuka_iiwa_arm/iiwa_common.h"
#include "drake/multibody/parsers/urdf_parser.h"
#include "drake/multibody/rigid_body_tree.h"

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

const int kDof = 7;
const int T = 50;
const int numStates = 14;
const double dt = 0.05;
const double finite_differences_epsilon = 1e-4;
const int maxIterations = 10;
const double converganceThreshold = 0.000001;
const double terminalPosWeight = 0.0001;
const double terminalVelWeight = 0.0001;

class iLQR {
 public:
  Eigen::MatrixXd k_;
  Eigen::MatrixXd K_;
  Eigen::VectorXd x0_;
  Eigen::MatrixXd X_;
  Eigen::VectorXd xt_;
  Eigen::MatrixXd U_;
  Eigen::MatrixXd Fm_;
  Eigen::MatrixXd fv_;
  Eigen::VectorXd optimal_cost_;
  Eigen::VectorXd current_cost_;
  double optimal_cost_sum_;
  double current_cost_sum_;

  bool update_rollout_ = false;
  double lambda_ = 1.0;
  double lambda_factor_ = 10.0;

  std::vector<MatrixXd> fx;
  std::vector<MatrixXd> fu;
  std::vector<MatrixXd> lxx;
  std::vector<MatrixXd> luu;
  std::vector<MatrixXd> lux;
  MatrixXd l;
  MatrixXd lx;
  MatrixXd lu;

  explicit iLQR(const RigidBodyTree<double>& tree)
      : tree_(tree)  {
    VerifyIiwaTree(tree);

    // traj_mat.resize(kNumJoints, T);
    l.resize(T,1);
    lx.resize(T,numStates);
    lu.resize(T,kDof);
  }

  void update(){
    optimal_cost_ = trajectoryCost(X_,U_);
    //TODO optimal_cost_sum_

    for(int i = 0; i < maxIterations; i++) {
      if(update_rollout_ == true){
        forwardPass(x0_,U_);
      }

      update_rollout_ = false;

      backwardsPass();
      // current_cost_ , Xc_ , Uc_ = generateAndEvalNewTrajectory(x0_,X_,U_)

      //TODO current_cost_sum_
      if (current_cost_sum_ < optimal_cost_sum_) {
        lambda_ /= lambda_factor_;
        update_rollout_ = true;
        //TODO what do I reuse / save?

        optimal_cost_sum_ = current_cost_sum_;
        optimal_cost_ = current_cost_;

        if ((fabs((current_cost_sum_ - optimal_cost_sum_))/fabs(optimal_cost_sum_)) < converganceThreshold) {
          break;
        }

      } else {
        lambda_ *= lambda_factor_;
      }
    }
  }


  void computeControlSequence(VectorXd x0, MatrixXd X, MatrixXd U) {
      MatrixXd newX(T,numStates);
      newX.row(0) = x0;

      double cost = 0.0;

      for (int i = 0; i < (T - 1) ; i++) {
        newX.row(i) = forwardDynamics(X.row(i) , U.row(i), i);
        cost += computeL(newX.row(i+1), U.row(i)) * dt;
      }

      // return X, cost TODO fix return / pass through
  }

  void forwardPass(VectorXd x0, MatrixXd U) {
    update_rollout_ = false;

    MatrixXd X(T, kDof);
    VectorXd cost(kDof);

    // control_sequence_rollout(x0 , U , X, cost)
    fx.clear();
    fu.clear();
    lxx.clear();
    luu.clear();
    lux.clear();

    for(int i = 0; i < T; i++) {
      fx.push_back(MatrixXd(numStates,numStates));
      fu.push_back(MatrixXd(numStates,numStates));
      lxx.push_back(MatrixXd(numStates,numStates));
      luu.push_back(MatrixXd(kDof,kDof));
      lux.push_back(MatrixXd(kDof,numStates));

    }

    for (int i = 0; i < T; i++) {

      l(i) = computeCost( X.row(i) , U.row(i), lx, lxx.at(i), lu, luu.at(i), lux.at(i), i);

      MatrixXd A(numStates,numStates);
      MatrixXd B(numStates,kDof);

      // TODO Why multiply bt dt
      l(i) *= dt;
      lx.row(i) *= dt;
      lxx.at(i) *= dt;
      lu.row(i) *= dt;
      luu.at(i) *= dt;
      lux.at(i) *= dt;

      // finiteDifferences(X.row(i), U.row(i), i, A, B);
      fx.at(i) = Eigen::MatrixXd::Identity(numStates,numStates) + A * dt;
      fu.at(i) = B * dt;

    }
      // #TODO ADD IN FINAL TIMESTEP COST.
  }

   VectorXd forwardDynamics(VectorXd x, MatrixXd u, int t) {
     MatrixXd M;
     MatrixXd C;
     MatrixXd G;
     double q[kDof];
     double qd[kDof];

     for(int i = 0; i < kDof; i++) {
       q[i] = x(i);
       qd[i] = x(i+kDof);
     }

     computeMCG(q, qd, M, C, G);


     return x;

   }

  double computeL(VectorXd x, VectorXd u) {
    return u.squaredNorm();
  }


  double computeCost(VectorXd x, VectorXd u, MatrixXd &lx, MatrixXd &lxx, MatrixXd &lu, MatrixXd &luu, MatrixXd &lux, int i ) {
    double l = u.squaredNorm();
    lx.row(i) = Eigen::VectorXd::Zero(numStates);
    lxx = Eigen::MatrixXd::Zero(numStates,numStates);
    lu.row(i)  = 2.0 * u;
    luu = 2.0 * Eigen::MatrixXd::Identity(kDof,kDof);
    lux = Eigen::MatrixXd::Zero(kDof, numStates);
    return l;
  }

  double computeFinalCost(VectorXd x, VectorXd u, MatrixXd &lx, MatrixXd &lxx) {
    VectorXd q = x.block<1,kDof>(0,0);
    VectorXd qd = x.block<1,kDof>(0,kDof);
    VectorXd dis = q - xt_;

    for (int i =0 ; i < dis.size(); i++) {
      dis(i) = dis(i) + 180.0;
      dis(i) = fmod((double)dis(i) ,360.0);
      dis(i) = dis(i) - 180.0;
    }

    double l = terminalPosWeight * dis.squaredNorm() + terminalVelWeight * qd.squaredNorm();



    return l;
  }
  //
  //     def cost_final(self, x):
  //         """ the final state cost function """
  //         num_states = x.shape[0]
  //         l_x = np.zeros((num_states))
  //         l_xx = np.zeros((num_states, num_states))
  //
  //         wp = 1e4 # terminal position cost weight
  //         wv = 1e4 # terminal velocity cost weight
  //
  //         xy = self.arm.x
  //         xy_err = np.array([xy[0] - self.target[0], xy[1] - self.target[1]])
  //         l = (wp * np.sum(xy_err**2) +
  //                 wv * np.sum(x[self.arm.DOF:self.arm.DOF*2]**2))
  //
  //         l_x[0:self.arm.DOF] = wp * self.dif_end(x[0:self.arm.DOF])
  //         l_x[self.arm.DOF:self.arm.DOF*2] = (2 *
  //                 wv * x[self.arm.DOF:self.arm.DOF*2])
  //
  //         eps = 1e-4 # finite difference epsilon
  //         # calculate second derivative with finite differences
  //         for k in range(self.arm.DOF):
  //             veps = np.zeros(self.arm.DOF)
  //             veps[k] = eps
  //             d1 = wp * self.dif_end(x[0:self.arm.DOF] + veps)
  //             d2 = wp * self.dif_end(x[0:self.arm.DOF] - veps)
  //             l_xx[0:self.arm.DOF, k] = ((d1-d2) / 2.0 / eps).flatten()
  //
  //         l_xx[self.arm.DOF:self.arm.DOF*2, self.arm.DOF:self.arm.DOF*2] = 2 * wv * np.eye(self.arm.DOF)
  //
  //         # Final cost only requires these three values
  // return l, l_x, l_xx

  void computeMCG(double * q_in, double * qd_in, MatrixXd &M, MatrixXd &C, MatrixXd &G) {
    double *qptr = &q_in[0];
    Eigen::Map<Eigen::VectorXd> q(qptr, kDof);
    double *qdptr = &qd_in[0];
    Eigen::Map<Eigen::VectorXd> qd(qdptr, kDof);
    KinematicsCache<double> cache = tree_.doKinematics(q, qd);
    const RigidBodyTree<double>::BodyToWrenchMap no_external_wrenches;


    /* Note that the 'dynamics bias term' \f$ C(q, v, f_\text{ext}) \f$ can be
    *computed by simply setting \f$ \dot{v} = 0\f$.
    * Note also that if only the gravitational terms contained in \f$ C(q, v,
    *f_\text{ext}) \f$ are required, one can set \a include_velocity_terms to
    *false.
    * Alternatively, one can pass in a KinematicsCache created with \f$ v = 0\f$
    *or without specifying the velocity vector.
    */
    C = tree_.dynamicsBiasTerm(cache, no_external_wrenches, true);
    G = tree_.dynamicsBiasTerm(cache, no_external_wrenches, false);
    M = tree_.massMatrix(cache);
  }


  void backwardsPass(){
    MatrixXd V = l.row(T);
    MatrixXd Vx = lx.row(T);
    MatrixXd Vxx = lxx.at(T);

    MatrixXd Qx(numStates , 1);
    MatrixXd Qu(kDof , 1);
    MatrixXd Qxx(numStates,numStates);
    MatrixXd Qux(kDof,numStates);
    MatrixXd Quu(kDof,kDof);
    MatrixXd Quu_inv;

    for (int i = 0; i < T - 1; i++) {
        Qx = lx.row(i) + fx.at(i).transpose() * Vx;
        Qu = lu.row(i) + fu.at(i).transpose() * Vx;
        Qxx = lxx.at(i) + fx.at(i).transpose() * (Vxx * fx.at(i));
        Qux = lux.at(i) + fu.at(i).transpose() * (Vxx * fx.at(i));
        Quu = luu.at(i) + fu.at(i).transpose() * (Vxx * fu.at(i));

        Eigen::JacobiSVD<MatrixXd> svd(Quu, Eigen::ComputeThinU | Eigen::ComputeThinV);
        MatrixXd U = svd.matrixU();
        MatrixXd V = svd.matrixV();
        MatrixXd S = svd.singularValues();
        //TODO Remove negative values

        S += lambda_ * MatrixXd::Ones(S.rows(),S.cols());
        S = S.asDiagonal().inverse();
        Quu_inv = U * (S * V.transpose());
        k_.row(i) = - Quu_inv * Qu;
        K_.row(i) = - Quu_inv * Qux;

        Vx = Qx - K_.row(i).transpose() * (Quu * k_.row(i));
        Vxx = Qxx - K_.row(i).transpose() * (Quu * K_.row(i));
    }
  }

 private:

   VectorXd trajectoryCost(MatrixXd X, MatrixXd U) {
     VectorXd cost(kDof);

     return cost;
   }


   lcm::LCM lcm_;
   const RigidBodyTree<double>& tree_;
};

int do_main(int argc, const char* argv[]) {

  auto tree = std::make_unique<RigidBodyTree<double>>();
  parsers::urdf::AddModelInstanceFromUrdfFileToWorld(
      GetDrakePath() + "/examples/kuka_iiwa_arm/urdf/iiwa14_estimated_params_fixed_gripper.urdf",
      multibody::joints::kFixed, tree.get());


  iLQR ilqr(*tree);

  double q_in[7] = {0,-0.05,0,0,0,0,0};
  double qd_in[7] = {0,0,0,0,0,0,0};
  MatrixXd M;
  MatrixXd C;
  MatrixXd G;
  // MatrixXd G_old;
  //
  // for (int i = 0; i < 10; i ++) {
  //   ilqr.computeMCG(q_in,qd_in,M,C,G);
  //   q_in[1] += 0.01;
  //   if(i>=1){
  //     std::cout << "G " << G(1) - G_old(1) << std::endl;
  //   }
  //   G_old = G;
  // }
  ilqr.computeMCG(q_in,qd_in,M,C,G);

  std::cout << "C " << C << std::endl;
  std::cout << "G " << G << std::endl;
  std::cout << "M " << M << std::endl;

  ilqr.update();
  // RobotPlanRunner runner(*tree);
  // runner.Run();

  return 0;
}

}  // namespace
}  // namespace kuka_iiwa_arm
}  // namespace examples
}  // namespace drake

int main(int argc, const char* argv[]) {
  return drake::examples::kuka_iiwa_arm::do_main(argc, argv);
}
