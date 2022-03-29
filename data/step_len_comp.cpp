#include <iostream>
#include <cmath>
#include <string>
#include <omp.h>

#include "../src/vmc.h"

int main(int argc, char **argv) {
    double step_len_vals[] = {0.05, 0.1, 0.5, 1.0};

    int N = 100;
    int d = 3;

    long mc_cycles = pow(2, 23);
    double equi_fraction = pow(2.0, -4.0);

    omp_set_num_threads(4);
    #pragma omp parallel for
    for (int sli = 0; sli < 4; sli++) {
        std::string filename, path;

        double step_length = step_len_vals[sli];

        System *ind_ho = new NonInteracting(N, d, 1);
        MonteCarlo *met = new Metropolis(step_length);

        Solver solver(ind_ho, met, mc_cycles, equi_fraction, (unsigned int) time(NULL));

        for (double alpha = 0.3; alpha <= 0.71; alpha += 0.04) {
            filename = "data_alpha" + std::to_string(alpha) + ".csv";
            path = "./part_b_sl_comp/sl" + std::to_string(step_length) + "/";

            solver.solve(alpha);
            solver.write_results(filename, path);
        }

        delete ind_ho;
        delete met;
    }

    return 0;
}
