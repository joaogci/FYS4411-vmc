#include <iostream>

#include "vmc.h"

int main() {
    System *ind_ho = new NonInteracting(100, 3, 1);
    MonteCarlo *met_sampler = new Metropolis(0.5);
    
    Solver solver(ind_ho, met_sampler, (long) 1e5, 0.3, 2022);

    solver.set_optimization_params(1e-4, 1, 1e-3, 100000, 0.3);
    double opt_alpha = solver.optimize_var_params(0.1);
    printf("opt_alpha %f \n", opt_alpha);

    solver.set_solve_params(1e6, 0.15);
    solver.solve(opt_alpha);

    return 0;
}

