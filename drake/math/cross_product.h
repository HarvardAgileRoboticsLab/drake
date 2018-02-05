#pragma once

#include <Eigen/Dense>

#include "drake/common/constants.h"
#include "drake/common/eigen_types.h"

namespace drake {
namespace math {
template <typename Derived>
drake::Matrix3<typename Derived::Scalar> VectorToSkewSymmetric(
    const Eigen::MatrixBase<Derived>& p) {
  EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Eigen::MatrixBase<Derived>,
                                           drake::kSpaceDimension);
  drake::Matrix3<typename Derived::Scalar> ret;
  ret.setZero();
  ret(1) = p(2);
  ret(2) = -p(1);
  ret(3) = -p(2);
  ret(5) = p(0);
  ret(6) = p(1);
  ret(7) = -p(0);
  
  //ret << 0.0, -p(2), p(1), p(2), 0.0, -p(0), -p(1), p(0), 0.0;
  return ret;
}

}  // namespace math
}  // namespace drake
