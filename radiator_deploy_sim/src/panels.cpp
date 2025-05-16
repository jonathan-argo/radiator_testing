#include "panels.h"
#include "gen.h"
#include "hinges.h"
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <cmath>


void panels::calcDistances(const Eigen::Matrix<double, 5, 1>& theta) {

    panels::r_1 = { 
        0.5 * panels::width1 * std::cos(theta[0]),
        0.5 * panels::width1 * std::sin(theta[0]), 
        0
    };
        
    panels::r_2 = { 
        panels::width1 * std::cos(theta[0]) + 0.5 * panels::width2 * std::cos(theta[1]),
        panels::width1 * std::sin(theta[0]) + 0.5 * panels::width2 * std::sin(theta[1]),
        0
    };
        
    panels::r_3 = {
        panels::width1 * std::cos(theta[0]) + panels::width2 * std::cos(theta[1]) + 0.5 * panels::width3 * std::cos(theta[2]),
        panels::width1 * std::sin(theta[0]) + panels::width2 * std::sin(theta[1]) + 0.5 * panels::width3 * std::sin(theta[2]),
        0
    };
        
    panels::r_4 = { 
        panels::width1 * std::cos(theta[0]) + panels::width2 * std::cos(theta[1]) + panels::width3 * std::cos(theta[2]) + 0.5 * panels::width4 * std::cos(theta[3]),
        panels::width1 * std::sin(theta[0]) + panels::width2 * std::sin(theta[1]) + panels::width3 * std::sin(theta[2]) + 0.5 * panels::width4 * std::sin(theta[3]),
        0
    };
        
    panels::r_5 = {
        panels::width1 * std::cos(theta[0]) + panels::width2 * std::cos(theta[1]) + panels::width3 * std::cos(theta[2]) + panels::width4 * std::cos(theta[3]) + 0.5 * panels::width5 * std::cos(theta[4]),
        panels::width1 * std::sin(theta[0]) + panels::width2 * std::sin(theta[1]) + panels::width3 * std::sin(theta[2]) + panels::width4 * std::sin(theta[3]) + 0.5 * panels::width5 * std::sin(theta[4]),
        0
    };

    panels::r_a = {
        0,
        0, 
        0
    };
        
    panels::r_b = { 
        panels::width1 * std::cos(theta[0]),
        panels::width1 * std::sin(theta[0]),
        0
    };
        
    panels::r_c = {
        panels::width1 * std::cos(theta[0]) + panels::width2 * std::cos(theta[1]),
        panels::width1 * std::sin(theta[0]) + panels::width2 * std::sin(theta[1]),
        0
    };
        
    panels::r_d = { 
        panels::width1 * std::cos(theta[0]) + panels::width2 * std::cos(theta[1]) + panels::width3 * std::cos(theta[2]),
        panels::width1 * std::sin(theta[0]) + panels::width2 * std::sin(theta[1]) + panels::width3 * std::sin(theta[2]),
        0
    };
        
    panels::r_e = {
        panels::width1 * std::cos(theta[0]) + panels::width2 * std::cos(theta[1]) + panels::width3 * std::cos(theta[2]) + panels::width4 * std::cos(theta[3]),
        panels::width1 * std::sin(theta[0]) + panels::width2 * std::sin(theta[1]) + panels::width3 * std::sin(theta[2]) + panels::width4 * std::sin(theta[3]),
        0
    };
}

void panels::calcMomInert() {
    panels::I_a = 0.3333 * panels::mass1 * panels::width1 * panels::width1;
    panels::I_b = 0.3333 * panels::mass2 * panels::width2 * panels::width2;
    panels::I_c = 0.3333 * panels::mass3 * panels::width3 * panels::width3;
    panels::I_d = 0.3333 * panels::mass4 * panels::width4 * panels::width4;
    panels::I_e = 0.3333 * panels::mass5 * panels::width5 * panels::width5;
}

panels::SystemMatrix panels::calcAccAndReac(const Eigen::Matrix<double, 10, 1>& theta_dtheta, panels::forceSumCoef& forceSumCoef) {

    Eigen::Matrix<double, 5, 1> theta = theta_dtheta.segment<5>(0);
    Eigen::Matrix<double, 5, 1> dtheta = theta_dtheta.segment<5>(5);

    Eigen::Matrix<double, 5, 1> I;
    I << panels::I_a, panels::I_b, panels::I_c, panels::I_d, panels::I_e;

    Eigen::Matrix<double, 5, 1> m;
    m << panels::mass1, panels::mass2, panels::mass3, panels::mass4, panels::mass5;

    Eigen::Matrix<double, 5, 1> k;
    k << hinges::k_a, hinges::k_b, hinges::k_c, hinges::k_d, hinges::k_e;

    Eigen::Matrix<double, 5, 1> w;
    w << panels::width1, panels::width2, panels::width3, panels::width4, panels::width5;


    // A for state variable:: ddtheta1, ddtheta2, ddtheta3, ddtheta4, ddtheta5, Re_1x, Re_1y, Re_2x, Re_2y, Re_3x, Re_3y, Re_4x, Re_4y, Re_5x, Re_5y
    Eigen::Matrix<double, 15, 15> A;
    A.setZero();

    Eigen::Matrix<double, 15, 1> b;
    b.setZero();

    // Sum Forces x
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j <= i; ++j) {
            A(i, j) = m(j) * forceSumCoef.accCoefX(i, j);
        }
        if (i < 4) {
            A(i, (2 * i + 7)) = 1;
        }
        A(i, (2 * i + 5)) = -1;
    }

    for (int i = 0; i < 5; ++i) {
        b(i) = m(i) * forceSumCoef.constTermX(i);
    }

    // Sum Forces y
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j <= i; ++j) {
            A(i + 5, j) = m(j) * forceSumCoef.accCoefY(i, j);
        }
        if (i < 4) {
            A(i + 5, (2 * i + 8)) = 1;
        }
        A(i + 5, (2 * i + 6)) = -1;
    }

    for (int i = 0; i < 5; ++i) {
        b(i + 5) = m(i) * forceSumCoef.constTermY(i);
    }

    // Sum Moments
    A(10, 0) = I(0); // this is good
    A(10, 8) = w(0) * std::cos(theta(0)); // this is good
    A(10, 7) = -w(0) * std::sin(theta(0)); // this is good
    b(10) = 4 * k(0) * (gen::pi - theta(0)) + 4 * k(1) * (0.5 * gen::pi - theta(0) + theta(1)); // this is good

    A(11, 1) = I(1); // this is good
    A(11, 10) = w(1) * std::cos(theta(1)); // this is good
    A(11, 9) = -w(1) * std::sin(theta(1)); // this is good
    b(11) = -4 * k(1) * (0.5 * gen::pi - theta(0) + theta(1)) - 4 * k(2) * (0.5 * gen::pi + theta(1) - theta(2)); // this is good

    A(12, 2) = I(2); // this is good
    A(12, 12) = w(2) * std::cos(theta(2)); // this is good
    A(12, 11) = -w(2) * std::sin(theta(2)); // this is good
    b(12) = 4 * k(2) * (0.5 * gen::pi + theta(1) - theta(2)) + 4 * k(3) * (0.5 * gen::pi - theta(2) + theta(3)); // this is good

    A(13, 3) = I(3); // this is good
    A(13, 14) = w(3) * std::cos(theta(3)); // this is good
    A(13, 13) = -w(3) * std::sin(theta(3)); // this is good
    b(13) = -4 * k(3) * (0.5 * gen::pi - theta(2) + theta(3)) - 4 * k(4) * (0.5 * gen::pi + theta(3) - theta(4)); // this is good

    A(14, 4) = I(4); // this is good
    b(14) = 4 * k(4) * (0.5 * gen::pi + theta(3) - theta(4)); 

    // Replace conditional friction with continuous friction model for all angles
    for (int i = 0; i < 5; ++i) {
        // Use tanh for smooth transition of friction around zero velocity
        b(10 + i) += -hinges::mu_friction * tanh(100 * dtheta(i));
    }

    for (int i = 0; i < 5; ++i) {
        b(10 + i) += -hinges::b_damp * dtheta(i);
    }

    // Mechanical stops to prevent over-rotation
    /*
    if (theta(0) > panels::theta_max1) {
        b(10) += -gen::k_stop * (theta(0) - panels::theta_max1);
    } 

    if ((theta(1) - theta(0)) < 0) {
        b(11) += gen::k_stop * (theta(1) - theta(0));
    }

    if ((theta(1) - theta(2)) < 0) {
        b(12) += gen::k_stop * (theta(1) - theta(2));
    }

    if ((theta(3) - theta(2)) < 0) {
        b(13) += gen::k_stop * (theta(3) - theta(2));
    }

    if ((theta(3) - theta(4)) < 0) {
        b(14) += gen::k_stop * (theta(3) - theta(4));
    }
    */

    gen::checkSVD(A);

    panels::SystemMatrix system;

    system.A = A;
    system.b = b;

    return system;
    
}

const Eigen::Matrix<double, 10, 1> panels::rk4(const Eigen::Matrix<double, 5, 4>& acceleration_buffer, const Eigen::Matrix<double, 10, 1>& theta_dtheta_n) {
    Eigen::Matrix<double, 10, 1> state_n_1;
    double dt = gen::time_step;
    
    Eigen::Matrix<double, 5, 1> theta = theta_dtheta_n.segment<5>(0);
    Eigen::Matrix<double, 5, 1> dtheta = theta_dtheta_n.segment<5>(5);

    Eigen::Matrix<double, 5, 1> k1, k2, k3, k4;
    Eigen::Matrix<double, 5, 1> l1, l2, l3, l4;

    for (int i = 0; i < 5; ++i) {
        k1(i) = dt * dtheta(i);
        l1(i) = dt * acceleration_buffer(i, 0);
        k2(i) = dt * (dtheta(i) + 0.5 * l1(i));
        l2(i) = dt * acceleration_buffer(i, 1);
        k3(i) = dt * (dtheta(i) + 0.5 * l2(i));
        l3(i) = dt * acceleration_buffer(i, 2);
        k4(i) = dt * (dtheta(i) + l3(i));
        l4(i) = dt * acceleration_buffer(i, 3);
    }

    for (int i = 0; i < 5; ++i) {
        state_n_1(i) = theta(i) + 0.16667 * (k1(i) + 2 * k2(i) + 2 * k3(i) + k4(i));
        state_n_1(i + 5) = dtheta(i) + 0.16667 * (l1(i) + 2 * l2(i) + 2 * l3(i) + l4(i));
    }

    return state_n_1;

}

const Eigen::Matrix<double, 10, 1> panels::semiImplicitEuler(const Eigen::Matrix<double, 5, 1>& ddtheta_n, const Eigen::Matrix<double, 10, 1>& theta_dtheta_n) {
    Eigen::Matrix<double, 5, 1> theta = theta_dtheta_n.segment<5>(0);
    Eigen::Matrix<double, 5, 1> dtheta = theta_dtheta_n.segment<5>(5);
    double dt = gen::time_step;

    Eigen::Matrix<double, 10, 1> theta_dtheta_n_1;

    for (int i = 0; i < 5; ++i) {
        theta_dtheta_n_1(i + 5) = dtheta(i) + dt * ddtheta_n(i);
        theta_dtheta_n_1(i) = theta(i) + dt * theta_dtheta_n_1(i + 5);
    }

    return theta_dtheta_n_1;
}

panels::forceSumCoef panels::calcAccCoef(const Eigen::Matrix<double, 10, 1>& state) {
    Eigen::Matrix<double, 5, 1> theta = state.segment<5>(0);
    Eigen::Matrix<double, 5, 1> dtheta = state.segment<5>(5);

    Eigen::Matrix<double, 5, 1> w;
    w << panels::width1, panels::width2, panels::width3, panels::width4, panels::width5;

    Eigen::Matrix<double, 5, 1> m;
    m << panels::mass1, panels::mass2, panels::mass3, panels::mass4, panels::mass5;

    panels::forceSumCoef forceSumCoef;
    forceSumCoef.accCoefX.setZero();
    forceSumCoef.accCoefY.setZero();
    forceSumCoef.constTermX.setZero();
    forceSumCoef.constTermY.setZero();

    double acc_coef_x = 0.5 * panels::width1 * (-std::sin(theta(0)));
    double curr_const_term_x = 0.5 * panels::width1 * dtheta(0) * dtheta(0) * std::cos(theta(0));

    double acc_coef_y = 0.5 * panels::width1 * (std::cos(theta(0)));
    double curr_const_term_y = 0.5 * panels::width1 * dtheta(0) * dtheta(0) * std::sin(theta(0));

    forceSumCoef.accCoefX(0, 0) = acc_coef_x;
    forceSumCoef.accCoefY(0, 0) = acc_coef_y;
    forceSumCoef.constTermX(0)  = curr_const_term_x;
    forceSumCoef.constTermY(0) = curr_const_term_y;

    double prev_const_x = curr_const_term_x;
    double prev_const_y = curr_const_term_y;

    for (int i = 1; i < 5; ++i) {
        forceSumCoef.accCoefX.row(i) = forceSumCoef.accCoefX.row(i - 1);
        forceSumCoef.accCoefX(i, i - 1) = 2 * forceSumCoef.accCoefX(i, i - 1);
        forceSumCoef.accCoefX(i, i) = 0.5 * w(i) * (-std::sin(theta(i))); // yeah these are correct

        forceSumCoef.accCoefY.row(i) = forceSumCoef.accCoefY.row(i - 1);
        forceSumCoef.accCoefY(i, i - 1) = 2 * forceSumCoef.accCoefY(i, i - 1);
        forceSumCoef.accCoefY(i, i) = 0.5 * w(i) * (std::cos(theta(i)));

        curr_const_term_x = 0.5 * w(i) * dtheta(i) * dtheta(i) * std::cos(theta(i));
        forceSumCoef.constTermX(i) = forceSumCoef.constTermX(i - 1) + prev_const_x + curr_const_term_x;
        prev_const_x = curr_const_term_x;

        curr_const_term_y = 0.5 * w(i) * dtheta(i) * dtheta(i) * std::sin(theta(i));
        forceSumCoef.constTermY(i) = forceSumCoef.constTermY(i - 1) + prev_const_y + curr_const_term_y;
        prev_const_y = curr_const_term_y;
    }

    return forceSumCoef;

}