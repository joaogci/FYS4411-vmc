#include <iostream>

#include "vmc.h"

int main() {
    System *ind_ho = new NonInteracting(100, 3, 1);
    MonteCarlo *met_sampler = new Metropolis(0.5);
    
    Solver solver(ind_ho, met_sampler, (long) 1e5, 0.3, 2022);

    solver.set_optimization_params(1e-4, 1, 1e-3, 100000, 0.3);
    printf("opt_alpha = %f \n", solver.optimize_var_params(0.1));

    // for (double alpha = 0.1; alpha <= 1.0; alpha+=0.1) {
    //     solver.solve(alpha);
    // }

    return 0;
}

