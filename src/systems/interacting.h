#ifndef NON_INTERACTING_H
#define NON_INTERACTING_H

#include "system.h"

class Interacting : public System {
private:

    long double distance(long double *rk, long double *rm) {
        return sqrt(SQUARE(rk[0] - rm[0]) + SQUARE(rk[1] - rm[1]) + SQUARE(rk[2] - rm[2]));
    }

public:

    Interacting(int N_, int dim_, double omega_) : System(N_, dim_, omega_) {}

    virtual long double evaluate_wf(long double *r) {
        double r2_sum = 0;
        for (int d = 0; d < dim; ++d) {
            r2_sum += SQUARE(r[d]);
        }
        
        return exp(- alpha * r2_sum);
    }

    virtual long double gradient_component_wf(long double x) {
        return 0;
    }

    virtual void gradient_wf(long double *grad, long double *rk) {
    }

    virtual long double laplacian_wf(long double *r) {
        double r2_sum = 0;
        for (int d = 0; d < dim; ++d) {
            r2_sum += SQUARE(r[d]);
        }

        return 4 * SQUARE(alpha) * r2_sum - 2 * alpha * dim;
    }

    virtual void quantum_force(long double *force, long double *r) {
        for (int d = 0; d < dim; d++) {
            force[d] = - 4 * alpha * r[d];
        }
    }

    virtual long double local_energy() {
        double r2_sum = 0;
        for (int i = 0; i < N; ++i) {
            for (int d = 0; d < dim; ++d) {
                r2_sum += SQUARE(r[i][d]);
            }
        }

        return N * dim * alpha - 2 * SQUARE(alpha) * r2_sum + 0.5 * SQUARE(omega) * r2_sum;
    }

};

#endif // NON_INTERACTING_H
