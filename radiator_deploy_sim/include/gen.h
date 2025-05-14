#ifndef GEN_H
#define GEN_H

#include <Eigen/Dense>
#include <iostream>

namespace gen {
    constexpr double pi = 3.14126535;
    constexpr double time_step = 0.01;
    constexpr double num_iter = 10000;

    void checkSVD(const Eigen::MatrixXd& A);
}

#endif // GEN_H