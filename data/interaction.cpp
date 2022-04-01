#include <iostream>
#include <cmath>

#include "../src/vmc.h"

int main(int argc, char **argv) {
    int N = atoi(argv[1]);
    int d = 3;
    double a = 0.00433;

    long mc_cycles = pow(2, 24);
    double equi_fraction = pow(2.0, -4.0);

    double eta = 1e-5;
    double tol = 1e-4;
    double h = 1e-3;
    long opt_cycles = pow(2, 20);
    double opt_fraction = pow(2.0, -2.0);

    System *ind_ho = new Interacting(N, d, 1, a);
    MonteCarlo *met = new Importance(0.1);

    Solver solver(ind_ho, met, mc_cycles, equi_fraction, (unsigned int) time(NULL));
    
    solver.set_optimization_params(eta, tol, h, opt_cycles, opt_fraction);
    double opt_alpha = solver.optimize_var_params(0.35);

    solver.solve(opt_alpha);

    solver.write_results("data_alpha" + std::to_string(opt_alpha) + "_N" + std::to_string(N) + ".csv", "./interaction/");

    delete ind_ho;
    delete met;

    return 0;
}
