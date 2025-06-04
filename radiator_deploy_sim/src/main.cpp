#include "panels.h"
#include "hinges.h"

#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <algorithm> // For std::clamp
#include <nlopt.hpp>

int main() {

    Eigen::Matrix<double, 2, 1> k;
    k << 
        0.003122977,
        0.003122977;

    Eigen::Matrix<double, 2, 1> theta = panels::simulate(k);

    std::cout << "System simulated with spring constants:\n" << k.transpose() << std::endl;

    return 0;
}