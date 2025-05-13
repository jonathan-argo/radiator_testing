#include "panels.h"
#include "gen.h"
#include "hinges.h"
#include <vector>
#include <iostream>
#include <cmath>


void panels::calcDistances(const std::vector<double>& theta) {

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

int panels::calcResidual(N_Vector y, N_Vector f, void *user_data) {
    double *y_data = N_VGetArrayPointer(y);
    double *f_data = N_VGetArrayPointer(f);

    double theta[5], dtheta[5], ddtheta[5], reaction_x[5], reaction_y[5];

    for (int i = 0; i < 5; ++i) {
        theta[i] = y_data[i];
        dtheta[i] = y_data[i + 5];
        ddtheta[i] = y_data[i + 10];
        reaction_x[i] = y_data[i + 15];
        reaction_y[i] = y_data[i + 20];
    };

    double I[5] = {panels::I_a, panels::I_b, panels::I_c, panels::I_d, panels::I_e};
    double m[5] = {panels::mass1, panels::mass2, panels::mass3, panels::mass4, panels::mass5};
    double k[5] = {hinges::k_a, hinges::k_b, hinges::k_c, hinges::k_d, hinges::k_e};
    double w[5] = {panels::width1, panels::width2, panels::width3, panels::width4, panels::width5};

    f_data[0] = (I[0] * ddtheta[0]) + 
                ((4 * k[0] - 4 * k[1]) * theta[0] + 4 * k[1] * theta[1]) + 
                (w[0] * std::cos(theta[0]) * reaction_y[1] - w[0] * std::sin(theta[0]) * reaction_x[1]) +
                (-4 * gen::pi * k[0] + 2 * gen::pi * k[1]);

    f_data[1] = (I[1] * ddtheta[1]) + 
                ((-4 * k[2] - 4 * k[1]) * theta[1] + 4 * k[1] * theta[0] + 4 * k[2] * theta[2]) + 
                (w[1] * std::cos(theta[1]) * reaction_y[2] - w[1] * std::sin(theta[1]) * reaction_x[2]) +
                (-2 * gen::pi * k[2] - 2 * gen::pi * k[1]);

    f_data[2] = (I[2] * ddtheta[2]) + 
                ((-4 * k[3] - 4 * k[2]) * theta[2] + 4 * k[2] * theta[1] - 4 * k[3] * theta[3]) + 
                (w[2] * std::cos(theta[2]) * reaction_y[3] - w[2] * std::sin(theta[2]) * reaction_x[3]) +
                (2 * gen::pi * k[2] - 2 * gen::pi * k[3]);

    f_data[3] = (I[3] * ddtheta[3]) + 
                ((-4 * k[4] - 4 * k[3]) * theta[3] + 4 * k[3] * theta[2] + 4 * k[4] * theta[4]) + 
                (w[3] * std::cos(theta[3]) * reaction_y[4] - w[3] * std::sin(theta[3]) * reaction_x[4]) +
                (-2 * gen::pi * k[3] - 2 * gen::pi * k[4]);

    f_data[4] = (I[4] * ddtheta[4]) + 
                (4 * k[4] * theta[3] - 4 * k[4] * theta[4]) + 
                (2 * gen::pi * k[4]);

    for (int i = 0; i < 4; ++i) {
        f_data[i + 5] = (-m[i] * ddtheta[i] * 0.5 * w[i] * std::sin(theta[i])) - 
                        (m[i] * theta[i] * theta[i] * 0.5 * w[i] * std::cos(theta[i])) +
                        (reaction_x[i + 1] - reaction_x[i]);
    }

    f_data[9] = (-m[4] * ddtheta[4] * 0.5 * w[4] * std::sin(theta[4])) - 
                (m[4] * theta[4] * theta[4] * 0.5 * w[4] * std::cos(theta[4])) +
                (-reaction_x[4]);

    for (int i = 0; i < 4; ++i) {
        f_data[i + 10] = (m[i] * ddtheta[i] * 0.5 * w[i] * std::cos(theta[i])) - 
                         (m[i] * theta[i] * theta[i] * 0.5 * w[i] * std::sin(theta[i])) +
                         (reaction_y[i + 1] - reaction_y[i]);
    }

    f_data[14] = (m[4] * ddtheta[4] * 0.5 * w[4] * std::cos(theta[4])) - 
                 (m[4] * theta[4] * theta[4] * 0.5 * w[4] * std::sin(theta[4])) +
                 (-reaction_y[4]);

    

    return 0;
}
