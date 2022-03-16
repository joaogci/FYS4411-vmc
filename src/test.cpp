#include <iostream>

#include "vmc.h"

int main() {
    System *ind_ho = new NonInteracting(100, 3, 1);
    MonteCarlo *met_sampler = new Metropolis(0.5);
    
    Solver solver(ind_ho, met_sampler, 1e5, 0.3, 2022);

    solver.set_optimization_params(1e-4, 1, 1e-3, pow(2.0, 14), pow(2.0, -2));
    double opt_alpha = solver.optimize_var_params(0.1);

    solver.set_solve_params(pow(2.0, 20), pow(2.0, -3));
    solver.solve(opt_alpha);

    solver.write_results("file.csv", "./");

    return 0;
}
