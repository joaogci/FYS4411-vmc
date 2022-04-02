#include <iostream>
#include <cmath>
#include <omp.h>

#include "../src/vmc.h"

#define NUM_THREADS 4

int main(int argc, char **argv) {
    int N = atoi(argv[1]);
    int d = 3;
    double a = 0.00433;

    long mc_cycles = pow(2, 24);
    double equi_fraction = pow(2.0, -4.0);

    double eta = 5e-5;
    double tol = 1e-6;
    double h = 1e-3;
    long opt_cycles = pow(2, 18);
    double opt_fraction = pow(2.0, -2.0);

    System *ind_ho = new Interacting(N, d, 1, a);
    MonteCarlo *met = new Importance(0.001);

    Solver solver(ind_ho, met, mc_cycles, equi_fraction, (unsigned int) time(NULL));
    
    solver.set_optimization_params(eta, tol, h, opt_cycles, opt_fraction);
    double opt_alpha = solver.optimize_var_params(0.35);

    delete ind_ho;
    delete met;

    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel 
    {
        System *ind_ho = new Interacting(N, d, 1, a);
        MonteCarlo *met = new Importance(0.001);

        Solver solver(ind_ho, met, mc_cycles / NUM_THREADS, equi_fraction, (unsigned int) time(NULL));
        solver.solve(opt_alpha);
        solver.write_results("data_interaction_thread" + std::to_string(omp_get_thread_num()) + "N" + std::string(N) + "_alpha" + std::to_string(opt_alpha) + ".csv", "./part_g_interaction/")
    
        delete ind_ho;
        delete met;
    }


    return 0;
}
