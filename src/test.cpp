#include <iostream>

#include "vmc.h"

int main() {
    System *ind_ho = new NonInteracting(100, 3, 1);
    MonteCarlo *met_sampler = new Metropolis(0.5);
    
    Solver solver(ind_ho, met_sampler, (long) 1e5, 0.3, 2022);

    for (double alpha = 0.1; alpha <= 1.0; alpha+=0.1) {
        solver.solve(alpha);
    }

    return 0;
}

