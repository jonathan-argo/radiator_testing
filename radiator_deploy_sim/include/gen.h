#ifndef GEN_H
#define GEN_H

#include <Eigen/Dense>
#include <iostream>

namespace gen {
    constexpr double pi = 3.14126535;

    void checkSVD(const Eigen::MatrixXd& A);
}

#endif // GEN_H