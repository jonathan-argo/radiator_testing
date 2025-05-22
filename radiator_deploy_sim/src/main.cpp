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

    Eigen::Matrix<double, 5, 1> k;
    k << 
        0.001561489,
        0.003122977,
        0.003122977,
        0.003122977,
        0.001561489;

    Eigen::Matrix<double, 5, 1> theta = panels::simulate(k);

    return 0;
}