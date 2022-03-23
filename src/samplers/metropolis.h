#ifndef METROPOLIS_H
#define METROPOLIS_H

#include "monte_carlo.h"

class Metropolis : public MonteCarlo {
private: 

    double step_length;

    virtual long double acceptence_ratio(long double *r_new, int k) {
        long double wf_new = system->evaluate_wf(r_new, k);
        long double wf_old = system->evaluate_wf(system->r[k], k);

        if (wf_old) {
            return SQUARE(wf_new) / SQUARE(wf_old);
        }
        return 0;
    }

    virtual void update_system(long double *r_new, int k) {
        for (int d = 0; d < system->dim; d++) {
            r_new[d] = system->r[k][d] + step_length * (rng->rand_uniform() - 0.5);
        }
    }

public:

    Metropolis(double step_length_) : MonteCarlo() {
        step_length = step_length_;
    }

};

#endif // METROPOLIS_H
