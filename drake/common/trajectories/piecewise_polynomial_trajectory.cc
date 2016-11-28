#include "drake/common/trajectories/piecewise_polynomial_trajectory.h"

#include <vector>

using Eigen::MatrixXd;

namespace drake {

PiecewisePolynomialTrajectory::PiecewisePolynomialTrajectory(
  const MatrixXd& trajectory_matrix, const std::vector<double>& times) {
  typedef PiecewisePolynomial<double> PPType;

  std::vector<PPType::PolynomialMatrix> polys;
  std::vector<double> segment_times;
  const int num_time_steps = times.size();
  // For each timestep, creates a PolynomialMatrix for each joint position.
  // Each column of trajectory_matrix represents a particular time, and the
  // rows of that column contain values for each joint coordinate.
  for (int i = 0; i < num_time_steps; ++i) {
    PPType::PolynomialMatrix poly_matrix(trajectory_matrix.rows(), 1);
    const auto traj_now = trajectory_matrix.col(i);

    Eigen::VectorXd vel_init(7);
    Eigen::VectorXd vel_final(7);

    // Produces cubic interpolating polynomials for each joint coordinate.
    if (i != num_time_steps - 1) {
      for (int row = 0; row < trajectory_matrix.rows(); ++row) {
        Eigen::Vector4d coeffs(0, 0, 0, 0);
        double t_inc_ = times.at(i+1) - times.at(i); // time duration of each trajectory segment

        // coarse estimation of desired velocities at the via points
        if (i == 1){
           vel_init(row) = 0;
           vel_final(row) = (trajectory_matrix(row, i + 2) - trajectory_matrix(row, i))/(times.at(i + 2) - times.at(i));
        }else if (i == num_time_steps - 3){
            vel_init(row) = (trajectory_matrix(row, i + 1) - trajectory_matrix(row, i - 1))/(times.at(i + 1) - times.at(i - 1));
            vel_final(row) = vel_init(row);
        }else{
           vel_init(row) = (trajectory_matrix(row, i + 1) - trajectory_matrix(row, i - 1))/(times.at(i + 1) - times.at(i - 1));
           vel_final(row) = (trajectory_matrix(row, i + 2) - trajectory_matrix(row, i))/(times.at(i + 2) - times.at(i));
        }
        // coefficients of cubic splines
        coeffs[0] = traj_now(row);
        coeffs[1] = vel_init(row); // coefficient of the first-order term
        coeffs[2] = 3.0*(trajectory_matrix(row, i + 1) - coeffs[0])/pow(t_inc_, 2) - 2.0*vel_init(row)/t_inc_ - vel_final(row)/t_inc_; // coefficient of the second-order term
        coeffs[3] = -2.0*(trajectory_matrix(row, i + 1) - coeffs[0])/pow(t_inc_, 3) + (vel_init(row) + vel_final(row))/pow(t_inc_, 2); // coefficient of the third-order term

        poly_matrix(row) = PPType::PolynomialType(coeffs);
      }
      polys.push_back(poly_matrix);
    }
    segment_times.push_back(times.at(i));
  }
  pp_ = PPType(polys, segment_times);
}

}  // namespace drake
