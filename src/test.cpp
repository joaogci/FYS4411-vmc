#include <iostream>

#include <omp.h>

#include "system.h"

int main(int argc, char **argv) {

    // System ind_ho((unsigned int) time(NULL));
    // System ind_ho(2);
    // ind_ho.init_particles(100, 3);
    // ind_ho.init_hamiltonian(1);
    // ind_ho.set_simulation_params(1000000L, 0.3, 0.5);

    // ind_ho.optimize_alpha(0.1, 1e-5, 1, 1e-4, 10000);   

    // for (double alpha = 0.1; alpha <= 1.0; alpha += 0.1) {
    //     ind_ho.run_metropolis(alpha);
    // }
    // ind_ho.write_results("res.csv", "./");

    long total_mc_cycles = 1e8;

    int n_cores = 4;
    omp_set_num_threads(n_cores);
    #pragma omp parallel
    {
        System ind_ho((unsigned int) (time(NULL)) * omp_get_thread_num());
        ind_ho.init_particles(100, 3);
        ind_ho.init_hamiltonian(1);
        ind_ho.set_simulation_params(total_mc_cycles / n_cores, 0.3, 0.5);
        ind_ho.run_metropolis(0.5);

        ind_ho.write_results("res" + std::to_string(omp_get_thread_num()) + ".csv", "./");        

    } 
    

    // printf("Total Wall Time: %fs \n", ind_ho.get_wall_time());

    return 0;
}
