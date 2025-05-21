#ifndef HINGES_H
#define HINGES_H

namespace hinges {
    inline double k_a;
    inline double k_b;
    inline double k_c;
    inline double k_d;
    inline double k_e;

    inline double prev_rad_force_a = 0;
    inline double prev_rad_force_b = 0;
    inline double prev_rad_force_c = 0;
    inline double prev_rad_force_d = 0;
    inline double prev_rad_force_e = 0;

    constexpr double b_damp = 0.01;
    constexpr double mu_friction = 0.1;
}

#endif // HINGES_H