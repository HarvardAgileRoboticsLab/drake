/// @file
///
/// Description!

#include <lcm/lcm-cpp.hpp>
#include <math.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include "drake/common/drake_assert.h"
#include "drake/common/drake_path.h"
#include "drake/examples/kuka_iiwa_arm/iiwa_common.h"
#include "drake/multibody/parsers/urdf_parser.h"
#include "drake/multibody/rigid_body_tree.h"
#include "drake/gps_run_controller.hpp"
#include "drake/gps_controller_gains.hpp"
#include "drake/lcmt_iiwa_status.hpp"
#include "bot_core/robot_state_t.hpp"
#include "robotlocomotion/robot_plan_t.hpp"
#include <boost/range/irange.hpp>
#include <chrono>
#include "drake/lcmt_iiwa_command.hpp"


using boost::irange;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using drake::Vector1d;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Dynamic;
using Eigen::RowMajor;
using Eigen::Matrix;
using Eigen::Ref;
using namespace std;
using ns = chrono::nanoseconds;
using get_time = chrono::steady_clock;



const char* const kLcmRunControllerChannel = "GPS_RUN_CONTROLLER";


namespace drake {
namespace examples {
namespace kuka_iiwa_arm {
namespace {

const int kDof = 7;
const int T = 200;
const int numStates = 14;
const double dt = 0.1;

const double finite_differences_epsilon = 1e-4;
const int maxIterations = 100;
const double converganceThreshold = 0.0;
const double terminalPosWeight = 1.0;
const double terminalVelWeight = 1.0;
const double min_energy_cost_weight = 1.0;
const double alpha = 1.0;
class iLQR {
 public:

  Eigen::Matrix<double,T,kDof,RowMajor> k_;
  std::vector< Matrix<double,kDof,numStates,RowMajor> > K_;
  Eigen::Matrix<double,numStates,RowMajor> x0_;
  Eigen::Matrix<double,T,numStates,RowMajor> X_;
  Eigen::Matrix<double,numStates,RowMajor> xt_;
  Eigen::Matrix<double,T,kDof,RowMajor> U_;
  double optimal_cost_;
  double current_cost_;


  bool update_rollout_ = true;
  double lambda_ = 1.0;
  double lambda_factor_ = 10.0;

  std::vector<Matrix<double,numStates,numStates,RowMajor>> fx_;
  std::vector<Matrix<double,numStates,kDof,RowMajor>> fu_;
  std::vector<Matrix<double,numStates,numStates,RowMajor>> lxx_;
  std::vector<Matrix<double,kDof,kDof,RowMajor>> luu_;
  std::vector<Matrix<double,kDof,numStates,RowMajor>> lux_;
  Matrix<double,T,RowMajor> l_;
  Matrix<double,T,numStates,RowMajor> lx_;
  Matrix<double,T,kDof,RowMajor> lu_;

  explicit iLQR(const RigidBodyTree<double>& tree)
      : tree_(tree)  {
    VerifyIiwaTree(tree);
    init_deratives();
    lcm_.subscribe(kLcmCommandChannel,
                    &iLQR::HandleCommand, this);
    lcm_.subscribe(kLcmStatusChannel,
                    &iLQR::HandleStatus, this);
  }

  void init_deratives(){
    l_ = Matrix<double,T,RowMajor>(T);
    lx_ = Matrix<double,T,numStates,RowMajor>(T,numStates);
    lu_ = Matrix<double,T,kDof,RowMajor>(T,kDof);;
    k_ = Matrix<double,T,kDof,RowMajor>(T,kDof);

    fx_.clear();
    fu_.clear();
    lxx_.clear();
    luu_.clear();
    lux_.clear();
    K_.clear();

    for(int i = 0; i < T; i++) {
      K_.push_back(Matrix<double,kDof,numStates,RowMajor>(kDof,numStates));
      fx_.push_back(Matrix<double,numStates,numStates,RowMajor>(numStates,numStates));
      fu_.push_back(Matrix<double,numStates,kDof,RowMajor>(numStates,kDof));
      lxx_.push_back(Matrix<double,numStates,numStates,RowMajor>(numStates,numStates));
      luu_.push_back(Matrix<double,kDof,kDof,RowMajor>(kDof,kDof));
      lux_.push_back(Matrix<double,kDof,numStates,RowMajor>(kDof,numStates));
    }
  }

  void update(){
    Eigen::Matrix<double,T,numStates,RowMajor> X_init = Eigen::Matrix<double,T,numStates,RowMajor>::Zero(T,numStates);
    computeControlSequence(X_.row(0) , U_, X_init);
    X_ = X_init;

    std::cout << "U " << U_  << ", X " << X_ <<std::endl;

    exit(0);
    optimal_cost_ = trajectoryCost(X_,U_);
    //TODO optimal_cost_sum_

    for(int i = 0; i < maxIterations; i++) {
      if(update_rollout_ == true){
        forwardPass(x0_,U_);
      }

      update_rollout_ = false;

      backwardsPass();
      // std::cout << "------------------------"<< std::endl;
      // std::cout << "------------------------" <<std::endl;
      // std::cout << "------------------------" <<std::endl;
      // std::cout << "------------------------" <<std::endl;

      Eigen::Matrix<double,T,kDof,RowMajor> newU = computeNewU(x0_);
      Eigen::Matrix<double,T,numStates,RowMajor> X_new = Eigen::Matrix<double,T,numStates,RowMajor>::Zero(T,numStates);
      current_cost_ = computeControlSequence(X_.row(0) , newU, X_new);

      std::cout << " i " << i << " Current Cost: " << current_cost_ << " Optimal Cost Sum: " <<  optimal_cost_ << std::endl;

      //TODO current_cost_sum_
      if (current_cost_ < optimal_cost_) {
        std::cout << "LOWER COST FOUND "  <<  std::endl;
        lambda_ /= lambda_factor_;
        update_rollout_ = true;
        //TODO what do I reuse / save?
        U_ = newU;
        optimal_cost_ = current_cost_;

        if ((fabs((current_cost_ - optimal_cost_))/fabs(optimal_cost_)) < converganceThreshold) {
          std::cout << "DONE!!!! " << std::endl;
          break;
        }

      } else {
        lambda_ *= lambda_factor_;
      }
    }
    Eigen::Matrix<double,T,numStates,RowMajor> X_new = Eigen::Matrix<double,T,numStates,RowMajor>::Zero(T,numStates);
    current_cost_ = computeControlSequence(X_.row(0) , U_, X_new);
    // std::cout << "k " << k_   <<std::endl;
    // std::cout << "k " << K_.at(T-1)   <<std::endl;

    std::cout << "U " << U_  << ", X " << X_new <<std::endl;
    std::cout << "update() METHOD DONE!!!! " << std::endl;
  }

  Matrix<double,T,kDof,RowMajor> computeNewU(Matrix<double,1,numStates,RowMajor> x0) {
    Matrix<double,T,kDof,RowMajor> U = Matrix<double,T,kDof,RowMajor>::Zero(T,kDof);
    Matrix<double,1,numStates,RowMajor> x = x0;

    for(int i =0; i < T;i++ ){
      U.row(i) = U_.row(i) + alpha * k_.row(i) + (K_.at(i) * (x - X_.row(i)).transpose()).transpose();
      // std::cout << U_.row(i) + alpha * k_.row(i) + (K_.at(i) * (x - X_.row(i)).transpose()).transpose() << std::endl;

      // std::cout << std::endl;
      // std::cout << "NEW U " << i * dt << " ,"<<  U.row(i) <<std::endl;
      // std::cout << "K " << K_.at(i) << std::endl;
      // std::cout << "k " << k_.row(i) << std::endl;
      // std::cout << "U " << U_.row(i) << std::endl;
      // // std::cout << "cost " << cost << std::endl;
      // std::cout << "x " << x << std::endl << std::endl;
      x = forwardDynamics(x,U.row(i));
      // std::cout << "x " << x << std::endl << std::endl;
      // std::cout << "X " << X_.row(i) << std::endl << std::endl;

      // std::cout << "1 " << U_.row(i) + k_.row(i) << std::endl << std::endl;
      // std::cout << "2 " << (x - X_.row(i)) << std::endl << std::endl;
      // std::cout << "3 " <<  K_.at(i) * (x - X_.row(i)) << std::endl << std::endl;
      //

    }
    return U;
  }

  double computeControlSequence(VectorXd x0, Matrix<double,T,kDof,RowMajor> U, Eigen::Ref<Matrix<double,T,numStates,RowMajor>> newX) {
      newX = Matrix<double,T,numStates,RowMajor>(T,numStates);
      newX.row(0).noalias() = x0;

      double cost = 0.0;

      for (int i = 0; i < (T - 1) ; i++) {
        newX.row(i+1).noalias() = forwardDynamics(newX.row(i) , U.row(i));
        cost += computeL(newX.row(i+1), U.row(i)) * dt;
        // std::cout << "newX " << i <<  newX.row(i)  << std::endl;
        // std::cout << "newX+1 " << newX.row(i+1)  << std::endl;
        // std::cout << "U " << U.row(i)  << std::endl;
        // std::cout << "c0: " << x0 << std::endl;
        // std::cout << "c0st: " << cost << std::endl;

        // std::cout << "-------------------------" << std::endl;
      }
      std::vector<Matrix<double,numStates,numStates,RowMajor>> tmplxx = lxx_;
      Matrix<double,T,numStates,RowMajor> tmplx = lx_;
      // std::cout << "-------------------------" << std::endl;
      // std::cout << "-------------------------" << newX.row(T-2)<<" ~~~~" << U.row(T-2) <<   std::endl;

      cost += computeFinalCost(newX.row(T-1), U.row(T-1), tmplx, tmplxx, T -1 ,terminalPosWeight,terminalVelWeight);
      // std::cout << "-------------------------" << std::endl;

      // std::cout << "-------------------------" << std::endl;
      // std::cout << "-------------------------" << std::endl;
      // std::cout << "-------------------------" << std::endl;
      // std::cout << "-------------------------" << std::endl;

      return cost;
  }

  void finiteDifferences(Ref<Matrix<double,numStates,numStates,RowMajor>> A, Ref<Matrix<double,numStates,kDof,RowMajor>> B, VectorXd x, VectorXd u, int t) {
      A = Matrix<double,numStates,numStates,RowMajor>(numStates,numStates);
      B = Matrix<double,numStates,kDof,RowMajor>(numStates, kDof);

      for(int i =0; i < numStates;i++) {
          Matrix<double,1,numStates,RowMajor> xPlus = x + Eigen::VectorXd::Ones(x.size()) *  finite_differences_epsilon;
          Matrix<double,1,numStates,RowMajor> xMinus = x -  Eigen::VectorXd::Ones(x.size()) *  finite_differences_epsilon;
          Matrix<double,1,numStates,RowMajor> xPlusNew = forwardDynamics(xPlus, u);
          // std::cout << "OLD X " << xPlus << " NEW X " << xPlusNew << " u " << u << std::endl;
          Matrix<double,1,numStates,RowMajor> xMinusNew = forwardDynamics(xMinus, u);
          Matrix<double,1,numStates,RowMajor> diff = (xPlusNew - xMinusNew) / (2 * finite_differences_epsilon);
          A.col(i).noalias() = diff;
      }
      for(int i =0; i < kDof;i++) {
          Matrix<double,1,kDof,RowMajor> uPlus = u + Eigen::VectorXd::Ones(u.size()) *  finite_differences_epsilon;
          Matrix<double,1,kDof,RowMajor> uMinus = u -  Eigen::VectorXd::Ones(u.size()) *  finite_differences_epsilon;
          Matrix<double,1,numStates,RowMajor> xPlusNew = forwardDynamics(x, uPlus);
          Matrix<double,1,numStates,RowMajor> xMinusNew = forwardDynamics(x, uMinus);
          Matrix<double,1,numStates,RowMajor> diff = (xPlusNew - xMinusNew) / (2 * finite_differences_epsilon);
          B.col(i).noalias() = diff;
      }
  }


  void forwardPass(VectorXd x0, Matrix<double,T,kDof,RowMajor> U) {
    update_rollout_ = false;

    // Matrix<double,T,kDof,RowMajor> X(T, kDof);
    Matrix<double,T,numStates,RowMajor> newX(T, numStates);
    // VectorXd cost(kDof);
    // std::cout << "forward pass A" << std::endl;
    computeControlSequence(x0 , U , newX);
    // std::cout << "forward pass B" << std::endl;

    X_ = newX;

    init_deratives();

    for (int i = 0; i < T-1; i++) {

      l_(i) = computeCost( newX.row(i) , U.row(i), lx_, lxx_, lu_, luu_, lux_, i);

      Matrix<double,numStates,numStates,RowMajor> A(numStates,numStates);
      Matrix<double,numStates,kDof,RowMajor> B(numStates,kDof);

      // TODO Why multiply bt dt
      // l_(i) *= dt;
      // lx_.row(i) *= dt;
      // lxx_.at(i) *= dt;
      // lu_.row(i) *= dt;
      // luu_.at(i) *= dt;
      // lux_.at(i) *= dt;

      finiteDifferences(A, B, newX.row(i), U.row(i), i);
      fx_.at(i).noalias() = Eigen::Matrix<double,numStates,numStates,RowMajor>::Identity(numStates,numStates) + A * dt;
      fu_.at(i).noalias() = B * dt;

      // std::cout << i << " lx" << lx_.row(i) << std::endl;
      // std::cout << i << " lxx" << lxx_.at(i) << std::endl;
      // std::cout << i << " lu" << lu_.row(i) << std::endl;
      // std::cout << i << " luu" << luu_.at(i) << std::endl;
      // std::cout << i << " lux" << lux_.at(i) << std::endl;
      // std::cout << i << " fx" << fx_.at(i) << std::endl;
      // std::cout << i << " fu" << fu_.at(i) << std::endl;

    }


    lu_.row(T-1).noalias()  = 0.0 * U.row(T-1);
    luu_.at(T-1).noalias() = 0.0 * Eigen::Matrix<double,kDof,kDof,RowMajor>::Identity(kDof,kDof);
    lux_.at(T-1).noalias() = Eigen::Matrix<double,kDof,numStates,RowMajor>::Identity(kDof, numStates);

    l_(T - 1)= computeFinalCost(newX.row(T-1), U.row(T-1), lx_, lxx_, T-1,terminalPosWeight,terminalVelWeight);
    // std::cout << "lu: " << l_(T-1) <<std::endl;
    // std::cout << "lxx" << lxx.at(T-1) << " , " << l(T-1)  << std::endl;
    fx_.at(T-1).noalias() = Eigen::Matrix<double,numStates,numStates,RowMajor>::Zero(numStates,numStates);
    fu_.at(T-1).noalias() = Eigen::Matrix<double,numStates,kDof,RowMajor>::Zero(numStates, kDof);

    // exit(0);
  }

   VectorXd forwardDynamics(VectorXd x, Matrix<double,1,kDof,RowMajor> u) {
     Matrix<double,Dynamic,Dynamic,RowMajor> M;
     Matrix<double,Dynamic,Dynamic,RowMajor> C;
     Matrix<double,Dynamic,Dynamic,RowMajor> G;
     VectorXd qd_0(kDof);
     VectorXd q_0(kDof);
     double q_in[kDof];
     double qd_in[kDof];

     for(int i = 0; i < kDof; i++) {
       q_in[i] = X_(T-2,i);
       qd_in[i] = X_(T-2,i+kDof);
       qd_0(i) = x(i+kDof);
       q_0(i) = x(i);
      //  std::cout << q_in[i] << " , " << q_0(i) << ": " << qd_in[i] << " , " << qd_0(i) << std::endl;
     }

     computeMCG(q_in, qd_in, M, C, G);

    //  std::cout << "---------------- C " << C << " G " << G.determinant() << " M " << M.determinant() << std::endl;

     Matrix<double,1,kDof,RowMajor> qdd;
     qdd.noalias() = M.inverse() * (u - C + G);
     VectorXd qd(kDof);
     qd.noalias() = qd_0 + qdd.transpose() * dt;
     VectorXd q(kDof);
     q.noalias() = q_0 + qd * dt + 0.5 * qdd.transpose() * pow(dt,2);
    //  std::cout << " Minv " << M.inverse() << " ,u " << u << " ,C " << C << std::endl;
    //  std::cout << " qd_0 " << qd_0 << " q_0 " << q_0 << "  x, " << x << std::endl;
    //  std::cout << " qdd " << qdd << " ,qd " << qd << " ,q " << q <<  "--------------" <<std::endl;

     VectorXd returnVec(2*kDof);
     returnVec << q,qd;
    //  std::cout << returnVec << std::endl;
     return returnVec;
   }

  double computeL(VectorXd x, VectorXd u) {
    return min_energy_cost_weight* u.squaredNorm();
  }


  double computeCost(VectorXd x, VectorXd u, Ref<Matrix<double,T,numStates,RowMajor>>lx, std::vector<Matrix<double,numStates,numStates,RowMajor>> &lxx, Ref<Matrix<double,T,kDof,RowMajor>> lu, std::vector<Matrix<double,kDof,kDof,RowMajor>> &luu, std::vector<Matrix<double,kDof,numStates,RowMajor>> &lux, int i ) {
    double l = min_energy_cost_weight * u.squaredNorm();
    lx.row(i).noalias() = Eigen::VectorXd::Zero(numStates);
    lxx.at(i).noalias() = Eigen::Matrix<double,numStates,numStates,RowMajor>::Zero(numStates,numStates);
    lu.row(i).noalias()  =  min_energy_cost_weight * 2.0 * u;
    luu.at(i).noalias() = 2.0 * min_energy_cost_weight * Eigen::Matrix<double,kDof,kDof,RowMajor>::Identity(kDof,kDof);
    lux.at(i).noalias() = Eigen::Matrix<double,kDof,numStates,RowMajor>::Zero(kDof, numStates);
    // return (l + computeFinalCost(x, u, lx, lxx, i ,terminalPosWeight/10.0,terminalPosWeight/10.0));
    return l;
  }

  double computeFinalCost(Matrix<double,1,numStates,RowMajor> x, Matrix<double,1,kDof,RowMajor> u, Ref<Matrix<double,T,numStates,RowMajor>> lx, std::vector<Matrix<double,numStates,numStates,RowMajor>> lxx, int t, double posWeight, double velWeight) {
    VectorXd q(7);// = x.block<1,kDof>(0,0);
    VectorXd qd(7);// = x.block<1,kDof>(0,kDof);
    VectorXd dis(7);// = q - xt_;

    for (int i =0 ; i < dis.size(); i++) { //TODO DOES THIS WORK IF I CHANGE IT TO RADS?
      q(i) = x(i);
      qd(i) = x(i+kDof);
      dis(i) = q(i) - xt_(i);
      dis(i) = dis(i) + (M_PI/ 2.0);
      dis(i) = fmod((double)dis(i) ,M_PI);
      dis(i) = dis(i) - (M_PI/ 2.0);
    }
    // std::cout << "DIStANCE " << dis<< " , x " << q << " , xt" << xt_ << std::endl;
    // lxx.at(t) = Eigen::Matrix<double,numStates,numStates,RowMajor>::Ones(numStates,numStates);
    double l = 0.5*  posWeight * dis.squaredNorm() + 0.5* velWeight * qd.squaredNorm();
    Eigen::Matrix<double,2*kDof,RowMajor> v(kDof * 2);
    // std::cout << "HERE" << std::endl ;

    for(int i =0; i < kDof; i++) {
      v(i) = posWeight * dis(i);
      v(kDof + i) = velWeight * qd(i);
    }
    // v << (2 * terminalPosWeight * dis) , (terminalVelWeight * qd);
    // std::cout << "DISTANCE " << dis << " Q " << q << " QD " << qd << " v " << v <<  std::endl;

    lx.row(t) = v;
    lxx.at(t) = posWeight * Eigen::Matrix<double,numStates,numStates,RowMajor>::Identity(numStates,numStates) + velWeight * Eigen::Matrix<double,numStates,numStates,RowMajor>::Identity(numStates,numStates);
    // std::cout << lxx.at(t) <<std::endl;
    return l;
  }


  void computeMCG(double * q_in, double * qd_in, Matrix<double,Dynamic,Dynamic,RowMajor> &M, Matrix<double,Dynamic,Dynamic,RowMajor> &C, Matrix<double,Dynamic,Dynamic,RowMajor> &G) {
    // for (int i = 0; i < kDof; i++){
    //   std::cout << " qin " << q_in[i] << " qd_in " << qd_in[i] << std::endl;
    //
    // }
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

    // C = Matrix<double, 1 , kDof,RowMajor>::Ones(1,kDof);
    // M = Matrix<double, kDof , kDof,RowMajor>::Identity(kDof,kDof);
    // G = Matrix<double, 1 , kDof,RowMajor>::Ones(1,kDof);


    // std::cout << "+++++++++++++++++++ C " << C << " G " << G << " M " << M << " Minv: " << M.inverse() << "------------------" << std::endl;
  }


  void backwardsPass(){
    // MatrixXd V = l_.row(T-1);

    Matrix<double,1,numStates,RowMajor> Vx;
    Vx.noalias() = lx_.row(T-1);
    Matrix<double,numStates,numStates,RowMajor> Vxx;
    Vxx.noalias() = lxx_.at(T-1);
    // std::cout << "Vxx" << Vxx.rows() << "," << Vxx.cols() << std::endl;
    // double V = l_(T-1);

    Matrix<double,1,numStates,RowMajor> Qx(1,numStates);
    Matrix<double,1,kDof,RowMajor> Qu(1,kDof);
    Matrix<double,numStates,numStates,RowMajor> Qxx(numStates,numStates);
    Matrix<double,kDof,numStates,RowMajor> Qux(kDof,numStates);
    Matrix<double,kDof,kDof,RowMajor> Quu(kDof,kDof);
    Matrix<double,kDof,kDof,RowMajor> Quu_inv(kDof,kDof);

    Matrix<double,kDof,kDof,RowMajor> Usvd(kDof,kDof);
    Matrix<double,kDof,kDof,RowMajor> Vsvd(kDof,kDof);
    Matrix<double,1,kDof,RowMajor> Ssvd(1,kDof);
    Eigen::JacobiSVD<Matrix<double,kDof,kDof,RowMajor>> svd;

    for (int i = T -1; i >= 0 ; i--) {
      // std::cout <<i << "," << lux_.size() << "," << Vxx.size() << " , " << fu_.size() << "," << fx_.size() << ":" << lux_.at(i) << std::endl << Vxx<< std::endl;

        Qx.noalias() = lx_.row(i) + (fx_.at(i).transpose() * Vx.transpose()).transpose();

        Qu.noalias() = lu_.row(i) + (fu_.at(i).transpose() * Vx.transpose()).transpose();
        Qxx.noalias() = lxx_.at(i) + fx_.at(i).transpose() * (Vxx * fx_.at(i));
        Qux.noalias() = lux_.at(i) + fu_.at(i).transpose() * (Vxx * fx_.at(i));
        Quu.noalias() = luu_.at(i) + fu_.at(i).transpose() * (Vxx * fu_.at(i));

        svd = Eigen::JacobiSVD<Matrix<double,kDof,kDof,RowMajor>>(Quu, Eigen::ComputeThinU | Eigen::ComputeThinV);
        Usvd.noalias() = svd.matrixU();
        Vsvd.noalias() = svd.matrixV();
        Ssvd.noalias() = svd.singularValues();

        // Eigen::EigenSolver<Matrix<double,kDof,kDof,RowMajor>> solver(Quu);
        // Matrix<double,Dynamic,Dynamic,RowMajor> values = solver.eigenvalues();
        // Matrix<double,Dynamic,Dynamic,RowMajor> vectors = solver.eigenvectors();

        for(int j = 0; j < Ssvd.rows(); j++) {
          for(int k = 0; k < Ssvd.cols();k++ ){

            if(Ssvd(j,k) < 0.0 ) {
              Ssvd(j,k) = 0.0;
            }else {
            }
            Ssvd(j,k) += lambda_;

            Ssvd(j,k) = 1.0 / Ssvd(j,k);

          }
        }

        // Quu_inv = vectors * (values.asDiagonal() * vectors.transpose());

        Quu_inv = Vsvd * (Ssvd.asDiagonal() * Usvd.transpose());

        k_.row(i) = - Quu_inv * Qu.transpose();
        K_.at(i) = - Quu_inv * Qux;

        Vx = Qx - (K_.at(i).transpose() * (Quu * k_.row(i).transpose())).transpose();
        Vxx = Qxx - K_.at(i).transpose() * (Quu * K_.at(i));


        if(i > T -5) {
        // // std::cout << "-------------" << i << "---------------" << std::endl;
        // std::cout << "Quu "  << Quu << std::endl;
        // // std::cout << "Qux "  << Qux << std::endl;
        // //
        // // if(i > 0) {
        // //   std::cout << "V: "  << (V) << std::endl;
        // //   V = V + l_(i-1);
        // // }
        // // std::cout << "lu "  << lu_.row(i) << std::endl;
        // // std::cout << "l "  << l_.row(i) << std::endl;
        // // std::cout << "lx "  << lx_.row(i) << std::endl;
        // // std::cout << "lxx "  << lxx_.at(i) << std::endl;
        // //
        // // std::cout << "delta_v "  << (-0.5 * k_.row(i) * Quu * k_.row(i).transpose()) << std::endl;
        // // std::cout << "V "  << Vsvd << std::endl;
        // // std::cout << "S "  << Ssvd << std::endl;
        // std::cout << "Quu_inv "  << Quu_inv << std::endl;
        // // std::cout << "k_ "  << k_.row(i) << std::endl;
        // // std::cout << "K_ "  << K_.at(i) << std::endl;
        // // std::cout << "Vx "  << Vx << std::endl;
        // std::cout << "Vxx "  << Vxx << std::endl;
        // std::cout << "Quu "  << Quu << std::endl;
        // std::cout << "Lux "  << lux_.at(i) << std::endl;
        //
        // std::cout << "fu  "  << fu_.at(i) << std::endl;
        // std::cout << "luu  "  << luu_.at(i) << std::endl;
        // std::cout << "fx  "  << fx_.at(i) << std::endl;

      }

    }
    // std::cout << "I AM DONE WITH THIS METHOD"<< std::endl;

  }


  void generate_init_data() {
     VectorXd cmd(7);
    //  send_reset(cmd);
    //  wait_for_convergance(cmd);
     VectorXd cmd2(7);
     for(int i = 0 ; i < 7; i++ ){
       cmd2(i) = 0.2;
     }

     send_reset(cmd2);

     chrono::steady_clock::time_point start = get_time::now();
     int currentStep = 0;

     while(!command_valid_ || !status_valid_) {
       lcm_.handle();
     }

     while (currentStep < T) {
       lcm_.handle();
       chrono::steady_clock::time_point end = get_time::now();
       chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(end - start);

       if (time_span.count() > dt) {
        //  std::cout << "dt: "  s<< time_span.count()  << std::endl;

         for(int i = 0; i < kDof; i++) {
           X_(currentStep,i) = iiwa_status_.joint_position_measured.at(i);
           if(currentStep != 0) {
             X_(currentStep,i+kDof) = (X_(currentStep,i) - X_((currentStep -1),i)) / time_span.count();
             if (X_(currentStep,i+kDof) == 0.0) {
               X_(currentStep,i+kDof) = X_(currentStep-1,i+kDof);
             }
           }
           U_(currentStep,i) = iiwa_command_.joint_torque.at(i);
          //  std::cout << iiwa_command_.joint_torque.at(i) << "," << iiwa_status_.joint_position_commanded.at(i) << std::endl;
         }

         currentStep++;
         start = get_time::now();
       }
     }

     std::cout << "X_ : " << X_ << std::endl;
     std::cout << "U_ : " << U_ << std::endl;

  }

 private:
   bool command_valid_ = false;
   bool status_valid_ = false;

   double trajectoryCost(Matrix<double,T,numStates,RowMajor> X, Matrix<double,T,kDof,RowMajor> U) {
     double cost = 0.0;

     for(int i =0; i < T -1 ;i++) {
       cost += computeL(X.row(i), U.row(i)) * dt;
     }
     std::vector<Matrix<double,numStates,numStates,RowMajor>> tmplxx = lxx_;
     Matrix<double,T,numStates,RowMajor> tmplx = lx_;
     cost += computeFinalCost(X.row(T-1), U.row(T-1), tmplx, tmplxx, T -1,terminalPosWeight,terminalVelWeight );
     return cost;
   }
   void HandleStatus(const lcm::ReceiveBuffer* rbuf, const std::string& chan,
                     const lcmt_iiwa_status* status) {
     iiwa_status_ = *status;
     status_valid_= true;
   }
   void HandleCommand(const lcm::ReceiveBuffer* rbuf, const std::string& chan,
                     const lcmt_iiwa_command* command) {
     iiwa_command_ = *command;
     command_valid_ = true;
   }


   void wait_for_convergance(VectorXd goal_pos){
     bool currently_resetting = true;
     double epsilon = 0.05;
     while (currently_resetting == true) {
       lcm_.handle();

       VectorXd curr(7);
       for (int i = 0; i < 7; i++) {
         curr(i) = iiwa_status_.joint_position_measured.at(i);
       }

       double diff = (goal_pos - curr).squaredNorm();

       if(diff < epsilon) {
         currently_resetting = false;
       }
     }
   }
   void send_reset(VectorXd goal_pos) {
       lcm_.handle();

       robotlocomotion::robot_plan_t plan;
       plan.utime = iiwa_status_.utime;
       plan.robot_name = "lbr_iiwa_14";
       plan.num_states = 2;
       plan.plan_info.push_back(1);
       plan.plan_info.push_back(1);

       bot_core::robot_state_t start;
       bot_core::robot_state_t goal;
       start.utime = 0.0;
       goal.utime = (T * dt * 1e6);
       start.num_joints = 7;
       goal.num_joints = 7;

       start.joint_name.push_back("iiwa_joint_1");
       start.joint_name.push_back("iiwa_joint_2");
       start.joint_name.push_back("iiwa_joint_3");
       start.joint_name.push_back("iiwa_joint_4");
       start.joint_name.push_back("iiwa_joint_5");
       start.joint_name.push_back("iiwa_joint_6");
       start.joint_name.push_back("iiwa_joint_7");
       goal.joint_name.push_back("iiwa_joint_1");
       goal.joint_name.push_back("iiwa_joint_2");
       goal.joint_name.push_back("iiwa_joint_3");
       goal.joint_name.push_back("iiwa_joint_4");
       goal.joint_name.push_back("iiwa_joint_5");
       goal.joint_name.push_back("iiwa_joint_6");
       goal.joint_name.push_back("iiwa_joint_7");

      for (int i = 0; i < 7; i++) {
        start.joint_velocity.push_back(0.0);
        goal.joint_velocity.push_back(0.0);
        start.joint_effort.push_back(0.0);
        goal.joint_effort.push_back(0.0);
        goal.joint_position.push_back(goal_pos(i));
        start.joint_position.push_back(iiwa_status_.joint_position_measured.at(i));
        // std::cout << "EEEE" << iiwa_status_.joint_position_measured.at(i) << "," << goal_pos(i)<< std::endl;
      }

      plan.plan.push_back(start);
      plan.plan.push_back(goal);

      lcm_.publish(kLcmPlanChannel, &plan);

   }
   lcmt_iiwa_command iiwa_command_;
   lcmt_iiwa_status iiwa_status_;
   lcm::LCM lcm_;
   const char* const kLcmStatusChannel = "IIWA_STATUS";
   const char* const kLcmPlanChannel = "COMMITTED_ROBOT_PLAN";
   const char* const kLcmCommandChannel = "IIWA_COMMAND";

   const RigidBodyTree<double>& tree_;
};

int do_main(int argc, const char* argv[]) {

  auto tree = std::make_unique<RigidBodyTree<double>>();
  parsers::urdf::AddModelInstanceFromUrdfFileToWorld(
      GetDrakePath() + "/examples/kuka_iiwa_arm/urdf/iiwa14_simplified_collision.urdf",
      multibody::joints::kFixed, tree.get());


  iLQR ilqr(*tree);
  ilqr.X_ = Eigen::Matrix<double,T,numStates,RowMajor>::Zero(T,numStates);
  ilqr.U_ = Eigen::Matrix<double,T,kDof,RowMajor>::Zero(T,kDof);
  ilqr.X_ += 0.0 * Eigen::Matrix<double,T,numStates,RowMajor>::Ones(T,numStates);
  ilqr.U_ += 0.0 * Eigen::Matrix<double,T,kDof,RowMajor>::Ones(T,kDof);

  ilqr.x0_ = Eigen::VectorXd(numStates) + 0.0 * Eigen::VectorXd::Ones(numStates);
  for(int i = kDof; i < 2* kDof; i ++) {
    ilqr.x0_(i) = 0.0;
  }
  ilqr.xt_.resize(14);

  for(int j = 0;j <kDof;j++){
    ilqr.xt_(j) = 0.25;
    ilqr.xt_(j+kDof) = 0.0;

  }
  ilqr.generate_init_data();
  ilqr.update();

  //
  // lcm::LCM lcm_;
  // gps_run_controller cmd;
  // cmd.numTimeSteps = 25;
  // cmd.numStates = 14;
  // cmd.dt = dt;
  //
  // for(int i =0; i <cmd.numTimeSteps; i++ ){
  //   gps_controller_gains K;
  //   gps_controller_gains k;
  //
  //   K.numStates = cmd.numStates;
  //   k.numStates = cmd.numStates;
  //
  //
  //   for(int j = 0; j < cmd.numStates; j++){
  //     // std::cout << "sdfsdfsdf"    << std::endl;
  //
  //     // K.values.push_back(ilqr.K_(i,j));
  //     k.values.push_back(ilqr.k_(i,j));
  //
  //   }
  //   cmd.K.push_back(K);
  //   cmd.k.push_back(k);
  //
  // }
  // lcm_.publish(kLcmRunControllerChannel, &cmd);
  //

  // double q_in[7] = {0,-0.05,0,0,0,0,0};
  // double qd_in[7] = {0,0,0,0,0,0,0};
  // Eigen::Matrix<double,Dynamic,Dynamic,RowMajor> M;
  // Eigen::Matrix<double,Dynamic,Dynamic,RowMajor> C;
  // Eigen::Matrix<double,Dynamic,Dynamic,RowMajor> G;
  // Eigen::Matrix<double,Dynamic,Dynamic,RowMajor> G_old;
  //
  // for (int i = 0; i < 10; i ++) {
  //   ilqr.computeMCG(q_in,qd_in,M,C,G);
  //   q_in[1] += 0.01;
  //   if(i>=1){
  //     std::cout << "G " << G(1) - G_old(1) << std::endl;
  //   }
  //   G_old = G;
  // }
  // ilqr.computeMCG(q_in,qd_in,M,C,G);

  // std::cout << "C " << C << std::endl;
  // std::cout << "G " << G << std::endl;
  // std::cout << "M " << M << std::endl;

  // ilqr.update();
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
