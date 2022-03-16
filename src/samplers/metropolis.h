#ifndef METROPOLIS_H
#define METROPOLIS_H

#include "monte_carlo.h"

class Metropolis : public MonteCarlo {
protected: 

    double step_length;

    virtual long double acceptence_ratio(long double *r_new, long double *r_old) {
        long double wf_new = system->evaluate_wf(r_new);
        long double wf_old = system->evaluate_wf(r_old);

        return SQUARE(wf_new) / SQUARE(wf_old);
    }

    virtual void update_system(long double *r_new, long double *r_old) {
        for (int d = 0; d < system->get_dim(); d++) {
            r_new[d] = r_old[d] + step_length * (rng->rand_uniform() - 0.5);
        }
    }

public:

    Metropolis(double step_length_) : MonteCarlo() {
        step_length = step_length_;
    }

};

#endif // METROPOLIS_H
