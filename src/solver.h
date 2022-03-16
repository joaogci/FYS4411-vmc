#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>

#include "rng.h"
#include "systems/system.h"
#include "samplers/monte_carlo.h"

class Solver {
private:

    RNG *rng;
    System* system;
    MonteCarlo* mc_sampler;

    long double *E;
    long double *E2;

    long mc_cycles;
    double equi_fraction;
    double measure_after;

public:

    Solver(System *system_, MonteCarlo* mc_sampler_, long mc_cycles_, double equi_fraction_, unsigned int seed_) {
        rng = new RNG(seed_);
        mc_sampler = mc_sampler_;
        system = system_;

        mc_sampler->set_rng(rng);
        mc_sampler->set_system(system);

        mc_cycles = mc_cycles_;
        equi_fraction = equi_fraction_;
        measure_after = mc_cycles * equi_fraction;
    }

    ~Solver() {
        delete rng;
        delete system;
        delete mc_sampler;
    }

    void solve(double alpha_) {
        int N = system->get_N();
        int dim = system->get_dim();
        int idx_p;
        long double r_new[dim];
        long double ratio;

        long double E_l, E = 0, E2 = 0;

        system->set_var_params(alpha_);

        for (int t = 0; t < mc_cycles; t++) {
            idx_p = rng->rand() % N;

            ratio = mc_sampler->step(r_new, system->get_rk(idx_p));
            if (ratio >= 1 || rng->rand_uniform() <= ratio) {
                system->update_rk(idx_p, r_new);
            }

            if(t >= measure_after) {
                E_l = system->local_energy();
                E += E_l;
                E2 += E_l * E_l;
            }
        }

        E /= (mc_cycles * (1.0 - equi_fraction));
        E2 /= (mc_cycles * (1.0 - equi_fraction));

        printf("Energy: %lf for alpha=%f \n", E, alpha_);
    }
    
    void optimize_var_params();

    void write_results();
    

};

#endif // SOLVER_H
