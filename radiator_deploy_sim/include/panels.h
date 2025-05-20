#ifndef PANELS_H
#define PANELS_H

#include "gen.h"
#include <vector> 
#include <Eigen/Dense>
#include <limits>
#include <algorithm> // For std::clamp
#include <fstream>

namespace panels {

// Variables (changeable)

    inline Eigen::Vector3d r_1;
    inline Eigen::Vector3d r_2;
    inline Eigen::Vector3d r_3;
    inline Eigen::Vector3d r_4;
    inline Eigen::Vector3d r_5;
    inline Eigen::Vector3d r_a;
    inline Eigen::Vector3d r_b;
    inline Eigen::Vector3d r_c; 
    inline Eigen::Vector3d r_d;
    inline Eigen::Vector3d r_e;
    inline Eigen::Vector3d r_tip;

    inline double I_a;
    inline double I_b;
    inline double I_c;
    inline double I_d;
    inline double I_e;

// Constants (unchangeable)

    // Panel 1 (root)
    constexpr double width1 = 0.385; // [m]
    constexpr double mass1 = 0.418; // [kg]
    constexpr double theta_init1 = 0.1; // [rad]
    constexpr double dtheta_init1 = 0.01; // [rad/s]
    constexpr double theta_max1 = gen::pi / 2; // [rad]

    // Panel 2
    constexpr double width2 = 0.720; // [m]
    constexpr double mass2 = 0.782; // [kg]
    constexpr double theta_init2 = gen::pi - 0.1; // [rad]
    constexpr double dtheta_init2 = -0.01; // [rad/s]

    // Panel 3
    constexpr double width3 = 0.720; // [m]
    constexpr double mass3 = 0.782; // [kg]
    constexpr double theta_init3 = 0.1; // [rad]
    constexpr double dtheta_init3 = 0.01; // [rad/s]

    // Panel 4
    constexpr double width4 = 0.720; // [m]
    constexpr double mass4 = 0.782; // [kg]
    constexpr double theta_init4 = gen::pi - 0.1; // [rad]
    constexpr double dtheta_init4 = -0.01; // [rad/s]

    // Panel 5
    constexpr double width5 = 0.720; // [m]
    constexpr double mass5 = 0.782; // [kg]
    constexpr double theta_init5 = 0.1; // [rad]
    constexpr double dtheta_init5 = 0.01; // [rad/s]

    const Eigen::Matrix<double, 5, 1> theta_init = (Eigen::Matrix<double, 5, 1>() << 
                theta_init1, theta_init2, theta_init3, theta_init4, theta_init5).finished();

    const Eigen::Matrix<double, 10, 1> init = (Eigen::Matrix<double, 10, 1>() << 
                theta_init1, theta_init2, theta_init3, theta_init4, theta_init5, 
                dtheta_init1, dtheta_init2, dtheta_init3, dtheta_init4, dtheta_init5).finished();



    struct SystemMatrix {
        Eigen::Matrix<double, 15, 15> A;
        Eigen::Matrix<double, 15, 1> b;
    };

    struct forceSumCoef {
        Eigen::Matrix<double, 5, 5> accCoefX;
        Eigen::Matrix<double, 5, 5> accCoefY;
        Eigen::Matrix<double, 5, 1> constTermX;
        Eigen::Matrix<double, 5, 1> constTermY;
    };

// Functions

    void calcDistances(const Eigen::Matrix<double, 5, 1>& theta);
    void calcMomInert();
    SystemMatrix calcAccAndReac(const Eigen::Matrix<double, 10, 1>& theta_dtheta, forceSumCoef& forceSumCoef);
    const Eigen::Matrix<double, 10, 1> semiImplicitEuler(const Eigen::Matrix<double, 5, 1>& ddtheta_n, const Eigen::Matrix<double, 10, 1>& theta_dtheta_n);
    const Eigen::Matrix<double, 10, 1> rk4(const Eigen::Matrix<double, 5, 4>& acceleration_buffer, const Eigen::Matrix<double, 10, 1>& theta_dtheta_n);
    forceSumCoef calcAccCoef(const Eigen::Matrix<double, 10, 1>& state);
    const Eigen::Matrix<double, 5, 1> simulate(const Eigen::Matrix<double, 5, 1>& k);
    
}

#endif // PANELS_H