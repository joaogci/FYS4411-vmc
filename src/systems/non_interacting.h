#ifndef NON_INTERACTING_H
#define NON_INTERACTING_H

#include "system.h"

class NonInteracting : public System{
public:

    NonInteracting(int N_, int dim_, double omega_) : System(N_, dim_, omega_) {}

    virtual void print_info() {
        printf("NonInteracting Bosons -> N: %d | d: %d | omega: %.3lf \n", N, dim, omega);
    }

    virtual long double evaluate_wf(long double *rk, int k) {
        long double tmp[dim];
        for (int d = 0; d < dim; d++) {
            tmp[d] = r[k][d];
            r[k][d] = rk[d];
        }

        long double r2_sum = 0;
        for (int i = 0; i < N; i++) {
            for (int d = 0; d < dim; ++d) {
                r2_sum += SQUARE(r[i][d]);
            }
        }

        for (int d = 0; d < dim; d++) {
            r[k][d] = tmp[d];
        }

        return exp(- alpha * r2_sum);
    }

    virtual void gradient_wf(long double *grad, long double *rk, int k) {
        for (int d = 0; d < dim; d++) {
            grad[d] = -2 * alpha * rk[d];
        }
    }

    virtual long double laplacian_wf(long double *rk, int k) {
        long double r2_sum = 0;
        for (int d = 0; d < dim; ++d) {
            r2_sum += SQUARE(rk[d]);
        }

        return 4 * SQUARE(alpha) * r2_sum - 2 * alpha * dim;
    }

    virtual void quantum_force(long double *force, long double *rk, int k) {
        for (int d = 0; d < dim; d++) {
            force[d] = - 4 * alpha * rk[d];
        }
    }

    virtual long double local_energy() {
        long double r2_sum = 0;
        for (int i = 0; i < N; ++i) {
            for (int d = 0; d < dim; ++d) {
                r2_sum += SQUARE(r[i][d]);
            }
        }

        return N * dim * alpha - 2 * SQUARE(alpha) * r2_sum + 0.5 * SQUARE(omega) * r2_sum;
    }

};

#endif // NON_INTERACTING_H
