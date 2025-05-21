#ifndef GEN_H
#define GEN_H

#include <Eigen/Dense>
#include <iostream>

namespace gen {
    constexpr double pi = 3.14159265359;
    constexpr double time_step = 0.01;
    constexpr double num_iter = 10000;
    constexpr double k_stop = 10000.0;
    const Eigen::VectorXd theta_target = (Eigen::VectorXd(5) << 0.5 * pi, 0.5 * pi, 0.5 * pi, 0.5 * pi, 0.5 * pi).finished(); // your target state
    const double alpha = 1e-3; // regularization weight

    void checkSVD(const Eigen::MatrixXd& A);
}

#endif // GEN_H