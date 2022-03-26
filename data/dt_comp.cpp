#include <iostream>
#include <cmath>
#include <string>
#include <omp.h>

#include "../src/vmc.h"

int main(int argc, char **argv) {
    double dt_vals[] = {10, 1, 0.1, 0.01, 0.001};

    int N = 100;
    int d = 3;

    long mc_cycles = pow(2, 23);
    double equi_fraction = pow(2.0, -4.0);

    omp_set_num_threads(4);
    #pragma omp parallel for
    for (int dti = 0; dti < 5; dti++) {
        std::string filename, path;

        double dt = dt_vals[dti];

        System *ind_ho = new NonInteracting(N, d, 1);
        MonteCarlo *met = new Importance(dt);

        Solver solver(ind_ho, met, mc_cycles, equi_fraction, (unsigned int) time(NULL));

        for (double alpha = 0.1; alpha < 1.0; alpha += 0.05) {
            filename = "data_alpha" + std::to_string(alpha) + "_dt" + std::to_string(dt) + ".csv";
            path = "./part_c/dt" + std::to_string(dt) + "/";

            solver.solve(alpha);
            solver.write_results(filename, path);
        }

        delete ind_ho;
        delete met;
    }

    return 0;
}
