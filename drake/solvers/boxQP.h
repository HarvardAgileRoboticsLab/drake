#pragma once

#include <Eigen/Dense>
#include <vector>

int boxQP(const Eigen::MatrixXd& Q, const Eigen::VectorXd& f,
        const Eigen::VectorXd& lb, const Eigen::VectorXd& ub,
        Eigen::VectorXd& x, Eigen::MatrixXd& Sfree, Eigen::MatrixXd& Pfree);
