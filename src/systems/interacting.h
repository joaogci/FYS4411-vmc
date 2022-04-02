#ifndef INTERACTING_H
#define INTERACTING_H

#include "system.h"

class Interacting : public System {
public:

    Interacting(int N_, int dim_, double omega_) : System(N_, dim_, omega_) {}

    Interacting(int N_, int dim_, double omega_, double a_) : System(N_, dim_, omega_, a_) {} 

    virtual void print_info() {
        printf("Interacting Bosons -> N: %d | d: %d | omega: %.3lf | omega_z: %.3f | a: %.6lf \n", N, dim, omega, beta, a);
    }

    virtual long double evaluate_wf(long double *rk, int k) {
        long double sing_part, jastrow;
        long double dist;
        double r2_sum = 0;

        long double tmp[dim];
        for (int d = 0; d < dim; d++) {
            tmp[d] = r[k][d];
            r[k][d] = rk[d];
        }

        for (int i = 0; i < N; i++) {
            for (int d = 0; d < dim; ++d) {
                r2_sum += (d != 2) ? SQUARE(r[i][d]) : beta * SQUARE(r[i][d]);
            }
        }
        sing_part = exp(- alpha * r2_sum);

        jastrow = 1.0;
        for (int i = 0; i < N && jastrow != 0; i++) {
            for (int j = 0; j < i && jastrow != 0; j++) {
                dist = DIST(r[j], r[i]);

                if (dist > a) {
                    jastrow *= 1.0 - a / dist;
                }
                else {
                    jastrow *= 0.0;
                }
            }
        }
        
        for (int d = 0; d < dim; d++) {
            r[k][d] = tmp[d];
        }

        return sing_part * jastrow;
    }

    virtual void gradient_wf(long double *grad, long double *rk, int k) {
        long double rkm;

        for (int d = 0; d < dim; ++d) {
            grad[d] = (d != 2) ? - 2 * alpha * rk[d] : - 2 * alpha * beta * rk[d];

            for (int m = 0; m < N; ++m) {
                if (m != k) {
                    rkm = DIST(rk, r[m]);
                    grad[d] += (rk[d] - r[m][d]) * (a / (SQUARE(rkm) * (rkm - a)));
                }
            }
        }
    }

    virtual long double laplacian_wf(long double *rk, int k) {
        long double line1 = - 2 * alpha * (dim - 1 + beta) + 
                            4 * SQUARE(alpha) * (SQUARE(rk[0]) + SQUARE(rk[1]) + SQUARE(beta) * SQUARE(rk[2]));
        long double sum1 = 0.0;
        long double sum2 = 0.0;
        long double sum3 = 0.0;

        long double rkm, rkn, denom1, denom2;
        long double rk_rm;

        for (int m = 0; m < N; ++m) {// Whitespace is a little bitch, @ me
            if (m != k) {
                rkm = DIST(rk, r[m]);
                denom1 = SQUARE(rkm) * (rkm - a);
                sum1 += (rk[0]*(rk[0] - r[m][0]) + rk[1]*(rk[1] - r[m][1]) + rk[2]*(rk[2] - r[m][2]))/(denom1);
                sum2 += ((dim-1)*a)/(denom1) + (a*a - 2*a*rkm)/(denom1 * (rkm - a));
            }
        }  

        long double line2 = - 4 * alpha * a * sum1;
        long double line3 = sum2;

        for (int m = 0; m < N; ++m) {
            for (int n = 0; n < N; ++n) {
                if (m != k && n != k) {
                    rk_rm = (rk[0] - r[m][0]) * (rk[0] - r[n][0])
                    + (rk[1] - r[m][1]) * (rk[1] - r[n][1])
                    + (rk[2] - r[m][2]) * (rk[2] - r[n][2]);

                    rkm = DIST(rk, r[m]);
                    rkn = DIST(rk, r[n]);

                    denom1 = SQUARE(rkm) * (rkm - a);
                    denom2 = SQUARE(rkn) * (rkn - a);
                    sum3 += rk_rm / (denom1 * denom2);
                }
            }
        }

        long double line4 = SQUARE(a) * sum3;

        return line1 + line2 + line3 + line4;
    }

    virtual void quantum_force(long double *force, long double *rk, int k) {
        gradient_wf(force, rk, k);
        for (int d = 0; d < dim; d++) {
            force[d] = 2 * force[d];
        }
    }

    virtual long double local_energy() {
        int k, d;
        long double kinetic_energy = 0;
        long double potential_energy = 0;

        for (k = 0; k < N; k++) {
            kinetic_energy += laplacian_wf(r[k], k);

            for (d = 0; d < dim; d++) {
                    potential_energy += (d != 2) ? SQUARE(r[k][d]) : SQUARE(beta) * SQUARE(r[k][d]);
            }
        }
        
        kinetic_energy = - kinetic_energy * 0.5;
        potential_energy = omega * potential_energy * 0.5;

        return potential_energy + kinetic_energy;
    }

};

#endif // INTERACTING_H
