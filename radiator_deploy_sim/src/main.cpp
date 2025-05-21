#include "panels.h"
#include "hinges.h"

#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <algorithm> // For std::clamp
#include <nlopt.hpp>

int main() {

    std::vector<double> k0{0.005, 0.005, 0.005, 0.005, 0.005};
    double minf;

    const int n = 5;
    nlopt::opt opt(nlopt::LN_NELDERMEAD, n);

    std::vector<double> lb(n, 0.001);
    std::vector<double> ub(n, 5);
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);

    opt.set_min_objective(panels::objective, nullptr);
    opt.set_xtol_rel(1e-4);
    opt.set_maxeval(2000);

    try {
        nlopt::result res = opt.optimize(k0, minf);
        std::cout << "Optimized stiffness values:\n";
        for (double ki : k0) std::cout << ki << " ";
        std::cout << "\nFinal cost: " << minf << "\n";
    } catch (std::exception& e) {
        std::cerr << "NLopt failed: " << e.what() << "\n";
    }

    return 0;
}