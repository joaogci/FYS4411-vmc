#include <iostream>
#include <cmath>
#include <string>
#include <omp.h>

#include "../src/vmc.h"

int main(int argc, char **argv) {
    int N_vals[4] = {1, 10, 100, 500};
    int d_vals[3] = {1, 2, 3};

    long mc_cycles = pow(2, 23);
    double equi_fraction = pow(2.0, -4.0);
    double step_length = 0.5;

    omp_set_num_threads(4);
    #pragma omp parallel for
    for (int Ni = 0; Ni < 4; Ni++) {
        #pragma omp parallel for 
        for (int di = 0; di < 3; di++) {
            std::string filename, path;

            int N = N_vals[Ni];
            int d = d_vals[di];

            System *ind_ho = new NonInteracting(N, d, 1);
            MonteCarlo *met = new Metropolis(step_length);

            Solver solver(ind_ho, met, mc_cycles, equi_fraction, (unsigned int) time(NULL));

            for (double alpha = 0.3; alpha <= 0.71; alpha += 0.04) {
                filename = "data_alpha" + std::to_string(alpha) + ".csv";
                path = "./part_b_test/N" + std::to_string(N) + "_d" + std::to_string(d) + "/";

                solver.solve(alpha);
                solver.write_results(filename, path);
            }

            delete ind_ho;
            delete met;
        }
    }

    return 0;
}
