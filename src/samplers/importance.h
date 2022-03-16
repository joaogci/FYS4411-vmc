#ifndef IMPORTANCE_H
#define IMPORTANCE_H

#include "monte_carlo.h"

class Importance : public MonteCarlo {
protected:

    double D = 0.5;
    double dt = 0.001;
    double sqrt_dt = sqrt(dt);

    virtual long double acceptence_ratio(long double *r_new, long double *r_old) {
        long double wf_new = system->evaluate_wf(r_new);
        long double wf_old = system->evaluate_wf(r_old);

        return greens_function(r_new, r_old) * SQUARE(wf_new) / SQUARE(wf_old);
    }

    virtual void update_system(long double *r_new, long double *r_old) {
        int dim = system->get_dim();
        long double q_force[dim];
        system->quantum_force(q_force, r_old);

        for (int d = 0; d < dim; d++) {
            r_new[d] = r_old[d] + D * q_force[d] * dt + sqrt_dt * rng->rand_normal();
        }
    }

    long double greens_function(long double *r_new, long double *r_old) {
        int dim = system->get_dim();
        long double gf = 0;
        long double q_force_old[dim];
        long double q_force_new[dim];
        system->quantum_force(q_force_old, r_old);
        system->quantum_force(q_force_new, r_new);

        for (int d = 0; d < dim; d++) {
            gf += 0.5 * (q_force_new[d] + q_force_old[d]) *
            (r_old[d] - r_new[d] + 0.5 * D * dt * (q_force_old[d] - q_force_new[d]));
        }

        return exp(gf);
    }

public:

    Importance() : MonteCarlo() {}

};

#endif // IMPORTANCE_H
