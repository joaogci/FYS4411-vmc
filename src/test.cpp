#include <iostream>

#include "system.h"

int main(int argc, char **argv) {

    System ind_ho((unsigned int) time(NULL));
    ind_ho.init_particles(100, 3);
    ind_ho.init_hamiltonian(1);
    ind_ho.set_simulation_params(1000000L, 0.2, 0.5);

    for (double alpha = 0.1; alpha <= 1.0; alpha += 0.1) {
        ind_ho.run_metropolis(alpha);
    }

    return 0;
}
