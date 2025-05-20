#include "panels.h"
#include "hinges.h"

#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <algorithm> // For std::clamp

int main() {

    Eigen::Matrix<double, 5, 1> k;
    k << 0.05, 0.04, 0.03, 0.02, 0.01;

    Eigen::Matrix<double, 5, 1> theta_10 = panels::simulate(k);

    std::cout << "State at 10 seconds:\n" << theta_10 << std::endl;
     
    return 0;
}