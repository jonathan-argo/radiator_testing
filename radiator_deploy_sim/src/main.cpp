#include "panels.h"
#include "hinges.h"

#include <Eigen/Dense>
#include <iostream>
#include <cmath>

int main() {

    hinges::k_a = 5;
    hinges::k_b = 4;
    hinges::k_c = 3;
    hinges::k_d = 2;
    hinges::k_e = 1;

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

    std::cout << "Root Moment of Inertia: " << panels::I_a << std::endl;
    std::cout << "Panel Moment of Intertia: " << panels::I_b << std::endl;

    panels::SystemMatrix system = panels::calcAccAndReac(panels::init);

    // sol in format: ddtheta1, ddtheta2, ddtheta3, ddtheta4, ddtheta5, Re_1x, Re_1y, Re_2x, Re_2y, Re_3x, Re_3y, Re_4x, Re_4y, Re_5x, Re_5y
    Eigen::VectorXd sol = system.A.fullPivLu().solve(system.b);

    std::cout << "Solution: \n" << sol << std::endl;

    return 0;
}