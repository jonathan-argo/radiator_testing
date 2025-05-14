#include "panels.h"
#include "hinges.h"

#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <cmath>

int main() {

    std::ofstream file("../state_sol.csv");
    if (!file.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        return 1;
    }

    file << "theta1,theta2,theta3,theta4,theta5,dtheta1,dtheta2,dtheta3,dtheta4,dtheta5,";
    file << "ddtheta1,ddtheta2,ddtheta3,ddtheta4,ddtheta5,Re_1x,Re_1y,Re_2x,Re_2y,Re_3x,Re_3y,Re_4x,Re_4y,Re_5x,Re_5y\n";

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

    panels::SystemMatrix system = panels::calcAccAndReac(state);

    // sol in format: ddtheta1, ddtheta2, ddtheta3, ddtheta4, ddtheta5, Re_1x, Re_1y, Re_2x, Re_2y, Re_3x, Re_3y, Re_4x, Re_4y, Re_5x, Re_5y
    Eigen::VectorXd sol = system.A.fullPivLu().solve(system.b);

    Eigen::Matrix<double, 5, 4> acceleration_buffer;

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

    Eigen::Matrix<double, 5, 1> ddtheta_n = sol.segment<5>(0);
    for (int col = 3; col > 0; --col) {
        acceleration_buffer.col(col) = acceleration_buffer.col(col - 1);
    }
    acceleration_buffer.col(0) = ddtheta_n;

    for (int i = 0; i < 3; ++i) {
        std::cout << "Acceleration Buffer: \n" << acceleration_buffer << std::endl;
        state = panels::semiImplicitEuler(ddtheta_n, state);
        panels::SystemMatrix system = panels::calcAccAndReac(state);
        Eigen::VectorXd sol = system.A.fullPivLu().solve(system.b);

        ddtheta_n << sol.segment<5>(0);
        for (int col = 3; col > 0; --col) {
            acceleration_buffer.col(col) = acceleration_buffer.col(col - 1);
        }
        acceleration_buffer.col(0) = ddtheta_n;

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

    std::cout << "Acceleration Buffer: \n" << acceleration_buffer << std::endl;

    double num_iter = gen::num_iter;
    for (int i = 0; i < num_iter - 3; ++i) {
        state = panels::rk4(acceleration_buffer, state);
        panels::SystemMatrix system = panels::calcAccAndReac(state);
        Eigen::VectorXd sol = system.A.fullPivLu().solve(system.b);

        ddtheta_n << sol.segment<5>(0);
        for (int col = 3; col > 0; --col) {
            acceleration_buffer.col(col) = acceleration_buffer.col(col - 1);
        }
        acceleration_buffer.col(0) = ddtheta_n;

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


    // System Matrices Debugging
    /*
    std::cout << "A: \n" << system.A << std::endl;
    std::cout << "b: \n" << system.b << std::endl;
    std::cout << "Solution: \n" << sol << std::endl;
    */

    

    return 0;
}