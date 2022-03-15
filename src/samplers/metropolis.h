#ifndef METROPOLIS_H
#define METROPOLIS_H

#include "monte_carlo.h"

class Metropolis : public MonteCarlo {
private: 

    double step_length;

public:

    Metropolis(System *system_, double step_length_, RNG *rng_) : MonteCarlo(system_, rng_) {
        step_length = step_length_;
    }

    long double acceptence_ratio(long double *r_new, long double *r_old) {
        long double wf_new = system->evaluate_wf(r_new);
        long double wf_old = system->evaluate_wf(r_old);

        return SQUARE(wf_new) / SQUARE(wf_old);
    }

    void update_system(long double *r_new, long double *r_old) {
        for (int d = 0; d < system->get_dim(); d++) {
            r_new[d] = r_old[d] + step_length * (rng->rand_uniform() - 0.5);
        }
    }

};

#endif // METROPOLIS_H
