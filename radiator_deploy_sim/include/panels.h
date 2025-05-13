#ifndef PANELS_H
#define PANELS_H

#include "gen.h"
#include <vector> 


namespace panels {

// Variables (changeable)

    inline std::vector<double> r_1;
    inline std::vector<double> r_2;
    inline std::vector<double> r_3;
    inline std::vector<double> r_4;
    inline std::vector<double> r_5;
    inline std::vector<double> r_a;
    inline std::vector<double> r_b;
    inline std::vector<double> r_c; 
    inline std::vector<double> r_d;
    inline std::vector<double> r_e;

    inline double I_a;
    inline double I_b;
    inline double I_c;
    inline double I_d;
    inline double I_e;

// Constants (unchangeable)

    // Panel 1 (root)
    constexpr double width1 = 0.385; // [m]
    constexpr double mass1 = 0.418; // [kg]
    constexpr double theta_init1 = 0; // [rad]

    // Panel 2
    constexpr double width2 = 0.720; // [m]
    constexpr double mass2 = 0.782; // [kg]
    constexpr double theta_init2 = gen::pi; // [rad]

    // Panel 3
    constexpr double width3 = 0.720; // [m]
    constexpr double mass3 = 0.782; // [kg]
    constexpr double theta_init3 = 0; // [rad]

    // Panel 4
    constexpr double width4 = 0.720; // [m]
    constexpr double mass4 = 0.782; // [kg]
    constexpr double theta_init4 = gen::pi; // [rad]

    // Panel 5
    constexpr double width5 = 0.720; // [m]
    constexpr double mass5 = 0.782; // [kg]
    constexpr double theta_init5 = 0; // [rad]

    static const std::vector<double> theta_init = {theta_init1, theta_init2, theta_init3, theta_init4, theta_init5};

// Functions

    void calcDistances(const std::vector<double>& theta);
    void calcMomInert();
    int calcResidual(N_Vector y, N_Vector f, void *user_data);
    
}

#endif // PANELS_H