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
    double measure_cycles;

    double eta = 1e-4;
    double tol = 1;
    double h = 1e-3;

    double solve_optimizer(double alpha_) {
        int N = system->get_N();
        int dim = system->get_dim();
        int idx_p;
        long double r_new[dim];
        long double ratio;

        long double E = 0;

        system->set_var_params(alpha_);

        for (int t = 0; t < mc_cycles; t++) {
            idx_p = rng->rand() % N;

            ratio = mc_sampler->step(r_new, system->get_rk(idx_p));
            if (ratio >= 1 || rng->rand_uniform() <= ratio) {
                system->update_rk(idx_p, r_new);
            }

            if(t >= measure_after) {
                E += system->local_energy();
            }
        }

        return E / measure_cycles;
    }

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
        measure_cycles = mc_cycles * (1.0 - equi_fraction);
    }

    ~Solver() {
        delete rng;
        delete system;
        delete mc_sampler;
    }

    void set_solve_params(long mc_cycles_, double equi_fraction_) {
        mc_cycles = mc_cycles_;
        equi_fraction = equi_fraction_;
        measure_after = mc_cycles * equi_fraction;
        measure_cycles = mc_cycles * (1.0 - equi_fraction);
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

        E /= measure_cycles;
        E2 /= measure_cycles;

        printf("Energy: %lf for alpha=%f \n", E, alpha_);
    }

    void set_optimization_params(double eta_, double tol_, double h_, long mc_cycles_, double equi_fraction_) {
        mc_cycles = mc_cycles_;
        equi_fraction = equi_fraction_;
        measure_after = mc_cycles * equi_fraction;
        measure_cycles = mc_cycles * (1.0 - equi_fraction); 

        eta = eta_;
        tol = tol_;
        h = h_;
    } 

    double optimize_var_params(double alpha_0_) {
        double opt_alpha = alpha_0_;
        long double derivative = tol + 1;
        long double Ep, Em;


        while (fabs(derivative) >= tol) {
            Ep = solve_optimizer(opt_alpha + h);
            Em = solve_optimizer(opt_alpha - h);

            derivative = (Ep - Em) / (2.0 * h);
            opt_alpha += - eta * derivative;
        }

        return opt_alpha;
    }

    void write_results();
    

};

#endif // SOLVER_H
