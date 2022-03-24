#include <iostream>
#include <cmath>
#include <string>
#include <omp.h>

#include "../src/vmc.h"

int main(int argc, char **argv) {
    int N_vals[4] = {1, 10, 100, 500};
    double step_len_vals[] = {0.05, 0.1, 0.5, 1.0};

    int d = 3;

    long mc_cycles = pow(2, 23);
    double equi_fraction = pow(2.0, -4.0);

    omp_set_num_threads(4);
    #pragma omp parallel for
    for (int Ni = 0; Ni < 4; Ni++) {
        for (int sli = 0; sli < 4; sli++) {
            std::string filename, path;

            int N = N_vals[Ni];
            double step_length = step_len_vals[sli];

            System *ind_ho = new NonInteracting(N, d, 1);
            MonteCarlo *met = new Metropolis(step_length);

            Solver solver(ind_ho, met, mc_cycles, equi_fraction, (unsigned int) time(NULL));

            for (double alpha = 0.1; alpha < 1.0; alpha += 0.05) {
                filename = "data_alpha" + std::to_string(alpha) + "_step_length" + std::to_string(step_length) + ".csv";
                path = "./non_interacting_step_len/N" + std::to_string(N) + "_d" + std::to_string(d) + "/";

                solver.solve(alpha);
                solver.write_results(filename, path);
            }

            delete ind_ho;
            delete met;
        }
    }

    return 0;
}
