#include <iostream>
#include <cmath>
#include <string>

#include "../src/vmc.h"

int main(int argc, char **argv) {
    int N = 100;
    int d = 3;

    long mc_cycles = pow(2, 23);
    double equi_fraction = pow(2.0, -13.0);
    double step_length = 0.5;

    std::string filename, path;

    System *ind_ho = new NonInteracting(N, d, 1);
    MonteCarlo *met = new Metropolis(step_length);

    Solver solver(ind_ho, met, mc_cycles, equi_fraction, (unsigned int) time(NULL));

    for (double alpha = 0.3; alpha <= 0.71; alpha += 0.04) {
         filename = "data_alpha" + std::to_string(alpha) + ".csv";
         path = "./part_b_test2/";

         solver.solve(alpha);
         solver.write_results(filename, path);
    }

    delete ind_ho;
    delete met;

    return 0;
}
