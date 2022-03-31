#include <iostream>
#include <cmath>

#include "../src/vmc.h"

int main(int argc, char **argv) {
    double alpha_0[5] = {0.2, 0.35, 0.5, 0.75, 0.8};

    int N = 100;
    int d = 3;

    long mc_cycles = pow(2, 22);
    double equi_fraction = pow(2.0, -2.0);

    double eta = 1e-5;
    double tol = 1e-4;
    double h = 1e-4;

    for (int a0 = 0; a0 < 5; a0++) {
        System *ind_ho = new NonInteracting(N, d, 1);
        MonteCarlo *met = new Metropolis(0.5);

        Solver solver(ind_ho, met, mc_cycles, equi_fraction, (unsigned int) time(NULL));

        solver.set_optimization_params(eta, tol, h, mc_cycles, equi_fraction);
        double opt_alpha = solver.optimize_var_params(alpha_0[a0]);

        delete ind_ho;
        delete met;
    }

    return 0;
}
