#include "panels.h"
#include "hinges.h"

#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <algorithm> // For std::clamp

// Function to check if a value is valid (not NaN or Inf)
bool isValidState(const Eigen::Matrix<double, 10, 1>& state) {
    for (int i = 0; i < state.size(); ++i) {
        if (std::isnan(state(i)) || std::isinf(state(i))) {
            return false;
        }
    }
    return true;
}

int main() {

    std::ofstream file("../data/state_sol.csv");
    if (!file.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        return 1;
    }

    std::ofstream file1("../data/distances.csv");
    if (!file1.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        return 1;
    }

    file << "time,theta1,theta2,theta3,theta4,theta5,dtheta1,dtheta2,dtheta3,dtheta4,dtheta5,";
    file << "ddtheta1,ddtheta2,ddtheta3,ddtheta4,ddtheta5,Re_1x,Re_1y,Re_2x,Re_2y,Re_3x,Re_3y,Re_4x,Re_4y,Re_5x,Re_5y\n";

    file1 << "time,r_ax,r_ay,r_bx,r_by,r_cx,r_cy,r_dx,r_dy,r_ex,r_ey,r_tipx,r_tipy\n";

    // Increased spring constants for better stability
    hinges::k_a = 0.05;
    hinges::k_b = 0.04;
    hinges::k_c = 0.03;
    hinges::k_d = 0.02;
    hinges::k_e = 0.01;

    panels::calcDistances(panels::theta_init);

    // Debugging Panel Distances
    /*
    std::cout << "Initial Distances: \n" << std::endl;
    std::cout << "r_1: \n" << panels::r_1 << std::endl; 
    std::cout << "r_2: \n" << panels::r_2 << std::endl; 
    std::cout << "r_3: \n" << panels::r_3 << std::endl; 
    std::cout << "r_4: \n" << panels::r_4 << std::endl; 
    std::cout << "r_5: \n" << panels::r_5 << std::endl; 
    std::cout << "r_a: \n" << panels::r_a << std::endl; 
    std::cout << "r_b: \n" << panels::r_b << std::endl; 
    std::cout << "r_c: \n" << panels::r_c << std::endl; 
    std::cout << "r_d: \n" << panels::r_d << std::endl; 
    std::cout << "r_e: \n" << panels::r_e << std::endl; 
    */

    panels::calcMomInert();

    // Moment of Inertia Debugging
    /*
    std::cout << "Root Moment of Inertia: " << panels::I_a << std::endl;
    std::cout << "Panel Moment of Intertia: " << panels::I_b << std::endl;
    */
    Eigen::Matrix<double, 10, 1> state = {
        panels::theta_init1,
        panels::theta_init2,
        panels::theta_init3,
        panels::theta_init4,
        panels::theta_init5,
        panels::dtheta_init1,
        panels::dtheta_init2,
        panels::dtheta_init3,
        panels::dtheta_init4,
        panels::dtheta_init5
    };

    Eigen::Matrix<double, 5, 1> theta = state.segment<5>(0);
    panels::calcDistances(theta);

    panels::forceSumCoef forceSumCoef;
    forceSumCoef = panels::calcAccCoef(state);

    //*
    std::cout << "Acceleration Coefficients X: \n" << forceSumCoef.accCoefX << std::endl;
    std::cout << "Acceleration Coefficients Y: \n" << forceSumCoef.accCoefY << std::endl;
    //*/

    panels::SystemMatrix system = panels::calcAccAndReac(state, forceSumCoef);

    // sol in format: ddtheta1, ddtheta2, ddtheta3, ddtheta4, ddtheta5, Re_1x, Re_1y, Re_2x, Re_2y, Re_3x, Re_3y, Re_4x, Re_4y, Re_5x, Re_5y
    Eigen::VectorXd sol = system.A.fullPivLu().solve(system.b);

    // System Matrices Debugging
    //*
    std::cout << "A: \n" << system.A << std::endl;
    std::cout << "b: \n" << system.b << std::endl;
    std::cout << "Solution: \n" << sol << std::endl;
    //*/

    Eigen::Matrix<double, 5, 4> acceleration_buffer;

    file << 0 << ",";
    file1 << 0 << ",";
    for (int j = 0; j < state.size(); ++j) {
            file << state(j) << ",";
        }

    // Write sol to CSV
    for (int j = 0; j < sol.size(); ++j) {
        file << sol(j);
        if (j < sol.size() - 1) {
            file << ",";
        }
    }

    file1 << panels::r_a.x() << ",";
    file1 << panels::r_a.y() << ",";
    file1 << panels::r_b.x() << ",";
    file1 << panels::r_b.y() << ",";
    file1 << panels::r_c.x() << ",";
    file1 << panels::r_c.y() << ",";
    file1 << panels::r_d.x() << ",";
    file1 << panels::r_d.y() << ",";
    file1 << panels::r_e.x() << ",";
    file1 << panels::r_e.y() << ",";
    file1 << panels::r_tip.x() << ",";
    file1 << panels::r_tip.y();
    file1 << "\n";

    Eigen::Matrix<double, 5, 1> ddtheta_n = sol.segment<5>(0);
    for (int col = 3; col > 0; --col) {
        acceleration_buffer.col(col) = acceleration_buffer.col(col - 1);
    }
    acceleration_buffer.col(0) = ddtheta_n;

    for (int i = 0; i < 3; ++i) {
        // std::cout << "Acceleration Buffer: \n" << acceleration_buffer << std::endl;
        state = panels::semiImplicitEuler(ddtheta_n, state);

        Eigen::Matrix<double, 5, 1> theta = state.segment<5>(0);
        panels::calcDistances(theta);

        file1 << gen::time_step * (1 + i) << ",";
        file1 << panels::r_a.x() << ",";
        file1 << panels::r_a.y() << ",";
        file1 << panels::r_b.x() << ",";
        file1 << panels::r_b.y() << ",";
        file1 << panels::r_c.x() << ",";
        file1 << panels::r_c.y() << ",";
        file1 << panels::r_d.x() << ",";
        file1 << panels::r_d.y() << ",";
        file1 << panels::r_e.x() << ",";
        file1 << panels::r_e.y() << ",";
        file1 << panels::r_tip.x() << ",";
        file1 << panels::r_tip.y();
        file1 << "\n";
        
        // Check for numerical instability in initial steps
        if (!isValidState(state)) {
            std::cerr << "Numerical instability detected in initial Euler step " << i << std::endl;
            std::cerr << "Current state: " << state.transpose() << std::endl;
            return 1; // Exit early if instability is detected
        }
        
        forceSumCoef = panels::calcAccCoef(state);
        panels::SystemMatrix system = panels::calcAccAndReac(state, forceSumCoef);
        Eigen::VectorXd sol = system.A.fullPivLu().solve(system.b);

        // Check solution quality
        double relative_error = (system.A * sol - system.b).norm() / system.b.norm();
        if (relative_error > 1e-6) {
            std::cerr << "Warning: Large relative error in initial system solution: " << relative_error << std::endl;
        }

        ddtheta_n << sol.segment<5>(0);
        for (int col = 3; col > 0; --col) {
            acceleration_buffer.col(col) = acceleration_buffer.col(col - 1);
        }
        acceleration_buffer.col(0) = ddtheta_n;

        file << gen::time_step * (1 + i) << ",";
        for (int j = 0; j < state.size(); ++j) {
            file << state(j) << ",";
        }

        // Write sol to CSV
        for (int j = 0; j < sol.size(); ++j) {
            file << sol(j);
            if (j < sol.size() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }

    // std::cout << "Acceleration Buffer: \n" << acceleration_buffer << std::endl;

    double num_iter = gen::num_iter;
    for (int i = 0; i < num_iter - 3; ++i) {
        state = panels::rk4(acceleration_buffer, state);

        Eigen::Matrix<double, 5, 1> theta = state.segment<5>(0);
        panels::calcDistances(theta);

        file1 << gen::time_step * (4 + i) << ",";
        file1<< panels::r_a.x() << ",";
        file1<< panels::r_a.y() << ",";
        file1<< panels::r_b.x() << ",";
        file1<< panels::r_b.y() << ",";
        file1<< panels::r_c.x() << ",";
        file1<< panels::r_c.y() << ",";
        file1<< panels::r_d.x() << ",";
        file1<< panels::r_d.y() << ",";
        file1<< panels::r_e.x() << ",";
        file1<< panels::r_e.y() << ",";
        file1<< panels::r_tip.x() << ",";
        file1<< panels::r_tip.y();
        file1 << "\n";
        
        // Check for numerical instability
        if (!isValidState(state)) {
            std::cerr << "Numerical instability detected at iteration " << i << std::endl;
            std::cerr << "Current state: " << state.transpose() << std::endl;
            break;
        }
        
        // Add bounds checking on angles (optional safety)
        for (int j = 0; j < 5; ++j) {
            // If angles exceed reasonable bounds, clamp them
            if (std::abs(state(j)) > 10.0) {
                // std::cerr << "Warning: Clamping theta[" << j << "] from " << state(j) << " to bounds" << std::endl;
                state(j) = std::clamp(state(j), -10.0, 10.0);
            }
        }
        
        forceSumCoef = panels::calcAccCoef(state);
        panels::SystemMatrix system = panels::calcAccAndReac(state, forceSumCoef);
        Eigen::VectorXd sol = system.A.fullPivLu().solve(system.b);

        // Check if solution is valid
        double relative_error = (system.A * sol - system.b).norm() / system.b.norm();
        if (relative_error > 1e-6) {
            std::cerr << "Warning: Large relative error in system solution: " << relative_error << std::endl;
        }

        ddtheta_n << sol.segment<5>(0);
        for (int col = 3; col > 0; --col) {
            acceleration_buffer.col(col) = acceleration_buffer.col(col - 1);
        }
        acceleration_buffer.col(0) = ddtheta_n;

        file << gen::time_step * (4 + i) << ",";
        for (int j = 0; j < state.size(); ++j) {
            file << state(j) << ",";
        }

        // Write sol to CSV
        for (int j = 0; j < sol.size(); ++j) {
            file << sol(j);
            if (j < sol.size() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }

    file.close();

    return 0;
}