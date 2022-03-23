#ifndef IMPORTANCE_H
#define IMPORTANCE_H

#include "monte_carlo.h"
#include <iostream>

class Importance : public MonteCarlo {
private:

    double D = 0.5;
    double dt = 0.001;
    double sqrt_dt = sqrt(dt);

    virtual long double acceptence_ratio(long double *r_new, int k) {
        long double wf_new = system->evaluate_wf(r_new, k);
        long double wf_old = system->evaluate_wf(system->r[k], k);

        if (wf_old) {
            return greens_function(r_new, system->r[k], k) * SQUARE(wf_new) / SQUARE(wf_old);
        }
        return 0;
    }

    virtual void update_system(long double *r_new, int k) {
        long double q_force[system->dim];
        system->quantum_force(q_force, system->r[k], k);

        for (int d = 0; d < system->dim; d++) {
            r_new[d] = system->r[k][d] + D * q_force[d] * dt + sqrt_dt * rng->rand_normal();
        }
    }

    long double greens_function(long double *r_new, long double *r_old, int k) {
        long double gf = 0;
        long double q_force_old[system->dim];
        long double q_force_new[system->dim];
        system->quantum_force(q_force_old, r_old, k);
        system->quantum_force(q_force_new, r_new, k);

        for (int d = 0; d < system->dim; d++) {
            gf += 0.5 * (q_force_new[d] + q_force_old[d]) *
            (r_old[d] - r_new[d] + 0.5 * D * dt * (q_force_old[d] - q_force_new[d]));
        }

        return exp(gf);
    }

public:

    Importance() : MonteCarlo() {}

};

#endif // IMPORTANCE_H
