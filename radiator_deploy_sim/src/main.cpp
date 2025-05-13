#include "panels.h"
#include "hinges.h"

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


    std::cout << "KINSOL version: " << SUNDIALS_VERSION_MAJOR << "." 
              << SUNDIALS_VERSION_MINOR << "." 
              << SUNDIALS_VERSION_PATCH << std::endl;


    panels::calcMomInert();




    return 0;
}