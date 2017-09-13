#pragma once

#include <string>

#include "drake/common/constants.h"
#include "drake/common/eigen_types.h"
#include "drake/math/autodiff.h"
#include "drake/math/roll_pitch_yaw.h"
#include "drake/multibody/joints/drake_joint_impl.h"
#include "drake/util/drakeGeometryUtil.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
class RollPitchYawFloatingJoint
    : public DrakeJointImpl<RollPitchYawFloatingJoint> {
 public:
  // disable copy construction and assignment
  // RollPitchYawFloatingJoint(const RollPitchYawFloatingJoint&) = delete;
  // RollPitchYawFloatingJoint& operator=(const RollPitchYawFloatingJoint&) =
  // delete;

 public:
  RollPitchYawFloatingJoint(const std::string& name,
                            const Eigen::Isometry3d& transform_to_parent_body)
      : DrakeJointImpl(*this, name, transform_to_parent_body, 6, 6) {}

  virtual ~RollPitchYawFloatingJoint() {}

  template <typename DerivedQ>
  Eigen::Transform<typename DerivedQ::Scalar, 3, Eigen::Isometry>
  jointTransform(const Eigen::MatrixBase<DerivedQ>& q) const {
    Eigen::Transform<typename DerivedQ::Scalar, 3, Eigen::Isometry> ret;
    auto pos = q.template middleRows<drake::kSpaceDimension>(0);
    auto rpy =
        q.template middleRows<drake::kRpySize>(drake::kSpaceDimension);
    ret.linear() = drake::math::rpy2rotmat(rpy);
    ret.translation() = pos;
    ret.makeAffine();
    return ret;
  }

  template <typename DerivedQ, typename DerivedMS>
  void motionSubspace(
      const Eigen::MatrixBase<DerivedQ>& q,
      // TODO(#2274) Fix NOLINTNEXTLINE(runtime/references).
      Eigen::MatrixBase<DerivedMS>& motion_subspace,
      typename drake::math::Gradient<DerivedMS, Eigen::Dynamic>::type*
          dmotion_subspace = nullptr) const {
    typedef typename DerivedQ::Scalar Scalar;
    motion_subspace.resize(drake::kTwistSize, get_num_velocities());
    auto rpy =
        q.template middleRows<drake::kRpySize>(drake::kSpaceDimension);
    Eigen::Matrix<Scalar, drake::kSpaceDimension, drake::kRpySize> E;
    rpydot2angularvelMatrix(rpy, E);
    Eigen::Matrix<Scalar, 3, 3> R = drake::math::rpy2rotmat(rpy);
    motion_subspace.template block<3, 3>(0, 0).setZero();
    motion_subspace.template block<3, 3>(0, 3) = R.transpose() * E;
    motion_subspace.template block<3, 3>(3, 0) = R.transpose();
    motion_subspace.template block<3, 3>(3, 3).setZero();

    if (dmotion_subspace) {
      dmotion_subspace->resize(motion_subspace.size(), get_num_positions());

      Scalar roll = rpy(0);
      Scalar pitch = rpy(1);
      Scalar yaw = rpy(2);

      Scalar cr = cos(roll);
      Scalar sr = sin(roll);
      Scalar cp = cos(pitch);
      Scalar sp = sin(pitch);
      Scalar cy = cos(yaw);
      Scalar sy = sin(yaw);

      dmotion_subspace->setZero();
      (dmotion_subspace->transpose())(112) = sr * sy + cr * cy * sp;
      (dmotion_subspace->transpose())(113) = cr*sy - cy*sp*sr;
      (dmotion_subspace->transpose())(118) = cr*sy - cy*sp*sr;
      (dmotion_subspace->transpose())(119) = -cr*cy - sp*sr*sy;
      (dmotion_subspace->transpose())(124) = cp*cr;
      (dmotion_subspace->transpose())(125) = -cp*sr;
      (dmotion_subspace->transpose())(133) = -sr;
      (dmotion_subspace->transpose())(134) = -cr;
      (dmotion_subspace->transpose())(139) = cp*cr;
      (dmotion_subspace->transpose())(140) = -cp*sr;
      (dmotion_subspace->transpose())(147) = -cy*sp;
      (dmotion_subspace->transpose())(148) = cp*cy*sr;
      (dmotion_subspace->transpose())(149) = cp*cr*cy;
      (dmotion_subspace->transpose())(153) = -sp*sy;
      (dmotion_subspace->transpose())(154) = cp*sr*sy;
      (dmotion_subspace->transpose())(155) = cp*cr*sy;
      (dmotion_subspace->transpose())(159) = -cp;
      (dmotion_subspace->transpose())(160) = -sp*sr;
      (dmotion_subspace->transpose())(161) = -cr*sp;
      (dmotion_subspace->transpose())(174) = -cp;
      (dmotion_subspace->transpose())(175) = -sp*sr;
      (dmotion_subspace->transpose())(176) = -cr*sp;
      (dmotion_subspace->transpose())(183) = -cp*sy;
      (dmotion_subspace->transpose())(184) = -cr*cy - sp*sr*sy;
      (dmotion_subspace->transpose())(185) = cy*sr - cr*sp*sy;
      (dmotion_subspace->transpose())(189) = cp*cy;
      (dmotion_subspace->transpose())(190) = cy*sp*sr - cr*sy;
      (dmotion_subspace->transpose())(191) = sr*sy + cr*cy*sp;


      // dmotion_subspace->transpose() << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      //     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      //     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      //     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      //     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      //     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      //     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      //     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      //     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, sr * sy + cr * cy * sp,
      //     cr * sy - cy * sp * sr, 0.0, 0.0, 0.0, 0.0, -cy * sr + cr * sp * sy,
      //     -cr * cy - sp * sr * sy, 0.0, 0.0, 0.0, 0.0, cp * cr, -cp * sr, 0.0,
      //     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -sr, -cr, 0.0, 0.0, 0.0, 0.0, cp * cr,
      //     -cp * sr, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -cy * sp, cp * cy * sr,
      //     cp * cr * cy, 0.0, 0.0, 0.0, -sp * sy, cp * sr * sy, cp * cr * sy,
      //     0.0, 0.0, 0.0, -cp, -sp * sr, -cr * sp, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      //     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -cp, -sp * sr, -cr * sp, 0.0, 0.0, 0.0,
      //     0.0, 0.0, 0.0, -cp * sy, -cr * cy - sp * sr * sy,
      //     cy * sr - cr * sp * sy, 0.0, 0.0, 0.0, cp * cy,
      //     -cr * sy + cy * sp * sr, sr * sy + cr * cy * sp, 0.0, 0.0, 0.0, 0.0,
      //     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      //     0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    }
  }

  template <typename DerivedQ, typename DerivedV>
  void motionSubspaceDotTimesV(
      const Eigen::MatrixBase<DerivedQ>& q,
      const Eigen::MatrixBase<DerivedV>& v,
      Eigen::Matrix<typename DerivedQ::Scalar, 6, 1>&
          motion_subspace_dot_times_v,
      typename drake::math::Gradient<
          Eigen::Matrix<typename DerivedQ::Scalar, 6, 1>, Eigen::Dynamic>::type*
          dmotion_subspace_dot_times_vdq = nullptr,
      typename drake::math::Gradient<
          Eigen::Matrix<typename DerivedQ::Scalar, 6, 1>, Eigen::Dynamic>::type*
          dmotion_subspace_dot_times_vdv = nullptr) const {
    typedef typename DerivedQ::Scalar Scalar;
    motion_subspace_dot_times_v.resize(drake::kTwistSize, 1);
    auto rpy =
        q.template middleRows<drake::kRpySize>(drake::kSpaceDimension);
    Scalar roll = rpy(0);
    Scalar pitch = rpy(1);
    Scalar yaw = rpy(2);

    auto pd = v.template middleRows<drake::kSpaceDimension>(0);
    Scalar xd = pd(0);
    Scalar yd = pd(1);
    Scalar zd = pd(2);

    auto rpyd =
        v.template middleRows<drake::kRpySize>(drake::kSpaceDimension);
    Scalar rolld = rpyd(0);
    Scalar pitchd = rpyd(1);
    Scalar yawd = rpyd(2);

    Scalar cr = cos(roll);
    Scalar sr = sin(roll);
    Scalar cp = cos(pitch);
    Scalar sp = sin(pitch);
    Scalar cy = cos(yaw);
    Scalar sy = sin(yaw);

    motion_subspace_dot_times_v[0] = -pitchd * yawd * cp;
    motion_subspace_dot_times_v[1] = rolld * yawd * cp * cr - pitchd * yawd *
      sp * sr - pitchd * rolld * sr;
    motion_subspace_dot_times_v[2] = -pitchd * rolld * cr - pitchd * yawd *
      cr * sp - rolld * yawd * cp * sr;
    motion_subspace_dot_times_v[3] = yd * (yawd * cp * cy - pitchd * sp * sy) -
      xd * (pitchd * cy * sp + yawd * cp * sy) - pitchd * zd * cp;
    motion_subspace_dot_times_v[4] = zd * (rolld * cp * cr - pitchd * sp * sr) +
      xd * (rolld * (sr * sy + cr * cy * sp) -
        yawd * (cr * cy + sp * sr * sy) + pitchd * cp * cy * sr) -
      yd * (rolld * (cy * sr - cr * sp * sy) +
        yawd * (cr * sy - cy * sp * sr) - pitchd * cp * sr * sy);
    motion_subspace_dot_times_v[5] = xd * (rolld * (cr * sy - cy * sp * sr) +
      yawd * (cy * sr - cr * sp * sy) + pitchd * cp * cr * cy) -
      zd * (pitchd * cr * sp + rolld * cp * sr) +
      yd * (yawd * (sr * sy + cr * cy * sp) -
        rolld * (cr * cy + sp * sr * sy) + pitchd * cp * cr * sy);

    if (dmotion_subspace_dot_times_vdq) {
      dmotion_subspace_dot_times_vdq->resize(motion_subspace_dot_times_v.rows(),
                                             get_num_positions());
      dmotion_subspace_dot_times_vdq->setZero();

      (dmotion_subspace_dot_times_vdq->transpose())(19) = -pitchd * rolld * cr - pitchd * yawd * cr * sp - rolld * yawd * cp * sr;
      (dmotion_subspace_dot_times_vdq->transpose())(20) = pitchd * rolld * sr + pitchd * yawd * sp * sr - rolld * yawd * cp * cr;
      (dmotion_subspace_dot_times_vdq->transpose())(22) = xd * (rolld * (cr * sy - cy * sp * sr) + yawd * (cy * sr - cr * sp * sy) + pitchd * cp * cr * cy) -
                                                          zd * (pitchd * cr * sp + rolld * cp * sr) + yd * (-rolld * (cr * cy + sp * sr * sy) +
                                                          yawd * (sr * sy + cr * cy * sp) + pitchd * cp * cr * sy);
      (dmotion_subspace_dot_times_vdq->transpose())(23) = -zd * (rolld * cp * cr - pitchd * sp * sr) -
                                                          xd * (rolld * (sr * sy + cr * cy * sp) -
                                                                yawd * (cr * cy + sp * sr * sy) + pitchd * cp * cy * sr) +
                                                          yd * (rolld * (cy * sr - cr * sp * sy) +
                                                                yawd * (cr * sy - cy * sp * sr) - pitchd * cp * sr * sy);
      (dmotion_subspace_dot_times_vdq->transpose())(24) = pitchd * yawd * sp;
      (dmotion_subspace_dot_times_vdq->transpose())(25) = -pitchd * yawd * cp * sr - rolld * yawd * cr * sp;
      (dmotion_subspace_dot_times_vdq->transpose())(26) = rolld * yawd * sp * sr - pitchd * yawd * cp * cr;
      (dmotion_subspace_dot_times_vdq->transpose())(27) = -xd * (pitchd * cp * cy - yawd * sp * sy) -
                                                           yd * (pitchd * cp * sy + yawd * cy * sp) + pitchd * zd * sp;
      (dmotion_subspace_dot_times_vdq->transpose())(28) = -zd * (pitchd * cp * sr + rolld * cr * sp) -
                                                          xd * (-rolld * cp * cr * cy + pitchd * cy * sp * sr +
                                                                yawd * cp * sr * sy) +
                                                          yd * (rolld * cp * cr * sy + yawd * cp * cy * sr -
                                                                pitchd * sp * sr * sy);
      (dmotion_subspace_dot_times_vdq->transpose())(29) = -zd * (pitchd * cp * cr - rolld * sp * sr) -
                                                          xd * (pitchd * cr * cy * sp + rolld * cp * cy * sr +
                                                                yawd * cp * cr * sy) -
                                                          yd * (-yawd * cp * cr * cy + pitchd * cr * sp * sy +
                                                                rolld * cp * sr * sy);
      (dmotion_subspace_dot_times_vdq->transpose())(33) = -xd * (yawd * cp * cy - pitchd * sp * sy) -
                                                           yd * (pitchd * cy * sp + yawd * cp * sy);
      (dmotion_subspace_dot_times_vdq->transpose())(34) = yd * (rolld * (sr * sy + cr * cy * sp) -
                                                          yawd * (cr * cy + sp * sr * sy) + pitchd * cp * cy * sr) +
                                                        xd * (rolld * (cy * sr - cr * sp * sy) +
                                                              yawd * (cr * sy - cy * sp * sr) - pitchd * cp * sr * sy);
      (dmotion_subspace_dot_times_vdq->transpose())(35) = yd * (rolld * (cr * sy - cy * sp * sr) +
                                                          yawd * (cy * sr - cr * sp * sy) + pitchd * cp * cr * cy) -
                                                        xd * (-rolld * (cr * cy + sp * sr * sy) +
                                                              yawd * (sr * sy + cr * cy * sp) + pitchd * cp * cr * sy);


      // dmotion_subspace_dot_times_vdq->transpose() << 0.0, 0.0, 0.0, 0.0, 0.0,
      //     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      //     -pitchd * rolld * cr - pitchd * yawd * cr * sp -
      //         rolld * yawd * cp * sr,
      //     pitchd * rolld * sr + pitchd * yawd * sp * sr -
      //         rolld * yawd * cp * cr,
      //     0.0,
      //     xd * (rolld * (cr * sy - cy * sp * sr) +
      //           yawd * (cy * sr - cr * sp * sy) + pitchd * cp * cr * cy) -
      //         zd * (pitchd * cr * sp + rolld * cp * sr) +
      //         yd * (-rolld * (cr * cy + sp * sr * sy) +
      //               yawd * (sr * sy + cr * cy * sp) + pitchd * cp * cr * sy),
      //     -zd * (rolld * cp * cr - pitchd * sp * sr) -
      //         xd * (rolld * (sr * sy + cr * cy * sp) -
      //               yawd * (cr * cy + sp * sr * sy) + pitchd * cp * cy * sr) +
      //         yd * (rolld * (cy * sr - cr * sp * sy) +
      //               yawd * (cr * sy - cy * sp * sr) - pitchd * cp * sr * sy),
      //     pitchd * yawd * sp, -pitchd * yawd * cp * sr - rolld * yawd * cr * sp,
      //     rolld * yawd * sp * sr - pitchd * yawd * cp * cr,
      //     -xd * (pitchd * cp * cy - yawd * sp * sy) -
      //         yd * (pitchd * cp * sy + yawd * cy * sp) + pitchd * zd * sp,
      //     -zd * (pitchd * cp * sr + rolld * cr * sp) -
      //         xd * (-rolld * cp * cr * cy + pitchd * cy * sp * sr +
      //               yawd * cp * sr * sy) +
      //         yd * (rolld * cp * cr * sy + yawd * cp * cy * sr -
      //               pitchd * sp * sr * sy),
      //     -zd * (pitchd * cp * cr - rolld * sp * sr) -
      //         xd * (pitchd * cr * cy * sp + rolld * cp * cy * sr +
      //               yawd * cp * cr * sy) -
      //         yd * (-yawd * cp * cr * cy + pitchd * cr * sp * sy +
      //               rolld * cp * sr * sy),
      //     0.0, 0.0, 0.0, -xd * (yawd * cp * cy - pitchd * sp * sy) -
      //                        yd * (pitchd * cy * sp + yawd * cp * sy),
      //     yd * (rolld * (sr * sy + cr * cy * sp) -
      //           yawd * (cr * cy + sp * sr * sy) + pitchd * cp * cy * sr) +
      //         xd * (rolld * (cy * sr - cr * sp * sy) +
      //               yawd * (cr * sy - cy * sp * sr) - pitchd * cp * sr * sy),

      //     yd * (rolld * (cr * sy - cy * sp * sr) +
      //           yawd * (cy * sr - cr * sp * sy) + pitchd * cp * cr * cy) -
      //         xd * (-rolld * (cr * cy + sp * sr * sy) +
      //               yawd * (sr * sy + cr * cy * sp) + pitchd * cp * cr * sy);
    }

    if (dmotion_subspace_dot_times_vdv) {
      dmotion_subspace_dot_times_vdv->resize(motion_subspace_dot_times_v.rows(),
                                             get_num_velocities());

      dmotion_subspace_dot_times_vdv->setZero();

      (dmotion_subspace_dot_times_vdv->transpose())(3) = -pitchd * cy * sp - yawd * cp * sy;
      (dmotion_subspace_dot_times_vdv->transpose())(4) = rolld * (sr * sy + cr * cy * sp) - yawd * (cr * cy + sp * sr * sy) +
                                                         pitchd * cp * cy * sr;
      (dmotion_subspace_dot_times_vdv->transpose())(5) = rolld * (cr * sy - cy * sp * sr) + yawd * (cy * sr - cr * sp * sy) +
                                                         pitchd * cp * cr * cy;
      (dmotion_subspace_dot_times_vdv->transpose())(9) = yawd * cp * cy - pitchd * sp * sy;
      (dmotion_subspace_dot_times_vdv->transpose())(10) = -rolld * (cy * sr - cr * sp * sy) - yawd * (cr * sy - cy * sp * sr) + pitchd * cp * sr * sy;
      (dmotion_subspace_dot_times_vdv->transpose())(11) = -rolld * (cr * cy + sp * sr * sy) + yawd * (sr * sy + cr * cy * sp) + pitchd * cp * cr * sy;
      (dmotion_subspace_dot_times_vdv->transpose())(15) = -pitchd * cp;
      (dmotion_subspace_dot_times_vdv->transpose())(16) = rolld * cp * cr - pitchd * sp * sr;
      (dmotion_subspace_dot_times_vdv->transpose())(17) = -pitchd * cr * sp - rolld * cp * sr;
      (dmotion_subspace_dot_times_vdv->transpose())(19) = -pitchd * sr + yawd * cp * cr;
      (dmotion_subspace_dot_times_vdv->transpose())(20) = -pitchd * cr - yawd * cp * sr;
      (dmotion_subspace_dot_times_vdv->transpose())(22) = xd * (sr * sy + cr * cy * sp) - yd * (cy * sr - cr * sp * sy) + zd * cp * cr;
      (dmotion_subspace_dot_times_vdv->transpose())(23) = xd * (cr * sy - cy * sp * sr) - yd * (cr * cy + sp * sr * sy) - zd * cp * sr;
      (dmotion_subspace_dot_times_vdv->transpose())(24) = -yawd * cp;
      (dmotion_subspace_dot_times_vdv->transpose())(25) = -sr * (rolld + yawd * sp);
      (dmotion_subspace_dot_times_vdv->transpose())(26) = -cr * (rolld + yawd * sp);
      (dmotion_subspace_dot_times_vdv->transpose())(27) = -zd * cp - xd * cy * sp - yd * sp * sy;
      (dmotion_subspace_dot_times_vdv->transpose())(28) = sr * (-zd * sp + xd * cp * cy + yd * cp * sy);
      (dmotion_subspace_dot_times_vdv->transpose())(29) = cr * (-zd * sp + xd * cp * cy + yd * cp * sy);
      (dmotion_subspace_dot_times_vdv->transpose())(30) = -pitchd * cp;
      (dmotion_subspace_dot_times_vdv->transpose())(31) = rolld * cp * cr - pitchd * sp * sr;
      (dmotion_subspace_dot_times_vdv->transpose())(32) = -pitchd * cr * sp - rolld * cp * sr;
      (dmotion_subspace_dot_times_vdv->transpose())(33) = cp * (yd * cy - xd * sy);
      (dmotion_subspace_dot_times_vdv->transpose())(34) = -xd * (cr * cy + sp * sr * sy) - yd * (cr * sy - cy * sp * sr);
      (dmotion_subspace_dot_times_vdv->transpose())(35) = xd * (cy * sr - cr * sp * sy) + yd * (sr * sy + cr * cy * sp);

      // dmotion_subspace_dot_times_vdv->transpose() << 0.0, 0.0, 0.0,
      //     -pitchd * cy * sp - yawd * cp * sy,
      //     rolld * (sr * sy + cr * cy * sp) - yawd * (cr * cy + sp * sr * sy) +
      //         pitchd * cp * cy * sr,
      //     rolld * (cr * sy - cy * sp * sr) + yawd * (cy * sr - cr * sp * sy) +
      //         pitchd * cp * cr * cy,
      //     0.0, 0.0, 0.0, yawd * cp * cy - pitchd * sp * sy,
      //     -rolld * (cy * sr - cr * sp * sy) - yawd * (cr * sy - cy * sp * sr) +
      //         pitchd * cp * sr * sy,
      //     -rolld * (cr * cy + sp * sr * sy) + yawd * (sr * sy + cr * cy * sp) +
      //         pitchd * cp * cr * sy,
      //     0.0, 0.0, 0.0, -pitchd * cp, rolld * cp * cr - pitchd * sp * sr,
      //     -pitchd * cr * sp - rolld * cp * sr, 0.0,
      //     -pitchd * sr + yawd * cp * cr, -pitchd * cr - yawd * cp * sr, 0.0,
      //     xd * (sr * sy + cr * cy * sp) - yd * (cy * sr - cr * sp * sy) +
      //         zd * cp * cr,
      //     xd * (cr * sy - cy * sp * sr) - yd * (cr * cy + sp * sr * sy) -
      //         zd * cp * sr,
      //     -yawd * cp, -sr * (rolld + yawd * sp), -cr * (rolld + yawd * sp),
      //     -zd * cp - xd * cy * sp - yd * sp * sy,
      //     sr * (-zd * sp + xd * cp * cy + yd * cp * sy),
      //     cr * (-zd * sp + xd * cp * cy + yd * cp * sy), -pitchd * cp,
      //     rolld * cp * cr - pitchd * sp * sr,
      //     -pitchd * cr * sp - rolld * cp * sr, cp * (yd * cy - xd * sy),
      //     -xd * (cr * cy + sp * sr * sy) - yd * (cr * sy - cy * sp * sr),
      //     xd * (cy * sr - cr * sp * sy) + yd * (sr * sy + cr * cy * sp);
    }
  }

  template <typename DerivedQ>
  void qdot2v(const Eigen::MatrixBase<DerivedQ>& q,
              Eigen::Matrix<typename DerivedQ::Scalar, Eigen::Dynamic,
                            Eigen::Dynamic, 0, DrakeJoint::MAX_NUM_VELOCITIES,
                        // TODO(#2274) Fix NOLINTNEXTLINE(runtime/references).
                            DrakeJoint::MAX_NUM_POSITIONS>& qdot_to_v,
              Eigen::Matrix<typename DerivedQ::Scalar, Eigen::Dynamic,
                            Eigen::Dynamic>* dqdot_to_v) const {
    qdot_to_v.setIdentity(get_num_velocities(), get_num_positions());
    drake::math::resizeDerivativesToMatchScalar(qdot_to_v, q(0));

    if (dqdot_to_v) {
      dqdot_to_v->setZero(qdot_to_v.size(), get_num_positions());
    }
  }

  template <typename DerivedQ>
  void v2qdot(const Eigen::MatrixBase<DerivedQ>& q,
              Eigen::Matrix<typename DerivedQ::Scalar, Eigen::Dynamic,
                            Eigen::Dynamic, 0, DrakeJoint::MAX_NUM_POSITIONS,
                        // TODO(#2274) Fix NOLINTNEXTLINE(runtime/references).
                            DrakeJoint::MAX_NUM_VELOCITIES>& v_to_qdot,
              Eigen::Matrix<typename DerivedQ::Scalar, Eigen::Dynamic,
                            Eigen::Dynamic>* dv_to_qdot) const {
    v_to_qdot.setIdentity(get_num_positions(), get_num_velocities());
    drake::math::resizeDerivativesToMatchScalar(v_to_qdot, q(0));

    if (dv_to_qdot) {
      dv_to_qdot->setZero(v_to_qdot.size(), get_num_positions());
    }
  }

  template <typename DerivedV>
  Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 1> frictionTorque(
      const Eigen::MatrixBase<DerivedV>& v) const {
    return Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 1>::Zero(
        get_num_velocities(), 1);
  }

  bool is_floating() const override { return true; }

// TODO(liang.fok) Remove this deprecated method prior to release 1.0.
#ifndef SWIG
  DRAKE_DEPRECATED("Please use is_floating().")
#endif
  bool isFloating() const override { return is_floating(); }

  Eigen::VectorXd zeroConfiguration() const override;
  Eigen::VectorXd randomConfiguration(
      std::default_random_engine& generator) const override;
  std::string get_position_name(int index) const override;

// TODO(liang.fok) Remove this deprecated method prior to release 1.0.
#ifndef SWIG
  DRAKE_DEPRECATED("Please use get_position_name().")
#endif
  std::string getPositionName(int index) const override;
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
#pragma GCC diagnostic pop  // pop -Wno-overloaded-virtual
