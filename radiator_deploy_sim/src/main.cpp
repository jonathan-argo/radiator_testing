#include "panels.h"
#include "hinges.h"

#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <algorithm> // For std::clamp
#include <nlopt.hpp>
#include <bitset>

int main() {

    constexpr int total_cases = 33; // 1 nominal + 32 boundaries
    int current_case = 0;

    // Progress bar setup (print before any simulation)
    auto update_progress = [](int current, int total) {
        float progress = float(current) / total;
        int bar_width = 50;
        std::cout << "\r[";
        int pos = int(bar_width * progress);
        for (int i = 0; i < bar_width; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0) << "% complete" << std::flush;
    };

    update_progress(current_case, total_cases);

    std::ofstream file("../data/state_sol.csv");
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open the file.");
    }

    file << "sim_id,time,theta1,theta2,theta3,theta4,theta5,dtheta1,dtheta2,dtheta3,dtheta4,dtheta5,";
    file << "ddtheta1,ddtheta2,ddtheta3,ddtheta4,ddtheta5,Re_1x,Re_1y,Re_2x,Re_2y,Re_3x,Re_3y,Re_4x,Re_4y,Re_5x,Re_5y\n";

    Eigen::Matrix<double, 5, 1> k;
    k << 
        0.0015,
        0.0015,
        0.0015,
        0.0015,
        0.0015;

    Eigen::Matrix<double, 5, 1> k_low;
    k_low = k * (1 - hinges::k_uncertainty);

    Eigen::Matrix<double, 5, 1> k_high;
    k_high = k * (1 + hinges::k_uncertainty);

    Eigen::Matrix<double, 5, 1> b;
    b << 
        hinges::b_damp,
        hinges::b_damp,
        hinges::b_damp,
        hinges::b_damp,
        hinges::b_damp;
        
    Eigen::Matrix<double, 5, 1> b_low;
    b_low = b * (1 - hinges::b_uncertainty);
    
    Eigen::Matrix<double, 5, 1> b_high;
    b_high = b * (1 + hinges::b_uncertainty);

    Eigen::Matrix<double, 5, 1> mu;
    mu <<
        hinges::mu_friction,
        hinges::mu_friction,
        hinges::mu_friction,
        hinges::mu_friction,
        hinges::mu_friction;
    
    Eigen::Matrix<double, 5, 1> mu_low;
    mu_low = mu * (1 - hinges::mu_uncertainty);

    Eigen::Matrix<double, 5, 1> mu_high;
    mu_high = mu * (1 + hinges::mu_uncertainty);

    std::bitset<6> sim_id("100000");

    Eigen::Matrix<double, 5, 1> theta = panels::simulate(k, sim_id, file);

    current_case++;
    update_progress(current_case, total_cases);

    for (int mask = 0; mask < (1 << 5); ++mask) {
        Eigen::Matrix<double, 5, 1> k_temp, mu_temp, b_temp;
        for (int i = 0; i < 5; ++i) {
            bool high = (mask >> i) & 1;
            k_temp[i] = high ? k_high[i] : k_low[i];
            mu_temp[i] = high ? mu_low[i] : mu_high[i];
            b_temp[i] = high ? b_low[i] : b_high[i];
        }

        std::bitset<6> label(mask);

        Eigen::Matrix<double, 5, 1> theta = panels::simulate(k_temp, label, file);
        
        current_case++;
        update_progress(current_case, total_cases);
    }

    file.close();
    std::cout << "\nSystem successfully simulated!" << std::endl;

    return 0;
}