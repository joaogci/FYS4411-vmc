#ifndef NON_INTERACTING_H
#define NON_INTERACTING_H

#include "system.h"

class Interacting : public System {
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
        long double rkm, denom;

        for (int d=0; d<dim; ++d){
            if (d<2){
                grad[d] = -2*alpha*rk[d];// + (d/(dim-1))*-2*alpha*beta*rk[d] <- Branchless programming attempt
            } else {
                grad[d] = -2*alpha*beta*rk[d];
            }
            for (int m=0; m<N; ++m){
                if (m != k){
                    rkm = DIST(rk, r[m]);
                    denom = SQUARE(rkm) * (rkm - a);
                    grad[d] += (rk[d] - r[m][d]) * (a/denom);
                }
            }
        }

    }

    virtual long double laplacian_wf(long double *rk, int k) {

        long double line1 = -2*alpha*(dim-1+beta) + 4*SQUARE(alpha) * (rk[0]*rk[0] + rk[1]*rk[1] + beta*rk[2]*rk[2]);
        long double sum1 = 0.0;
        long double sum2 = 0.0;
        long double sum3 = 0.0;

        long double rkm, rkn, rkm2, rkn2, denom1, denom2;

        for (int m=0; m<N; ++m)                                             {// Whitespace is a little bitch, @ me
            if (m != k){
                rkm = DIST(rk, r[m]);
                rkm2 = SQUARE(rkm);
                denom1 = rkm2*(rkm - a);
                sum1 += (rk[0]*(rk[0] - r[m][0]) + rk[1]*(rk[1] - r[m][1]) + rk[2]*(rk[2] - r[m][2]))/(denom1);
                sum2 += ((d-1)*a)/(denom1) + (a*a - 2*a*rkm)/(denom1 * (rkm - a));
            }
        }  

        long double line2 = -4*alpha*a *sum1;
        long double line3 = sum2;
        long double rk_rm;

        for (int m=0; m<N; ++m){
            for (int n=0; n<N; ++n) {
                if (m != k && n != k){
                    rk_rm = (rk[0] - r[m][0]) * (rk[0] - r[n][0])
                    + (rk[1] - r[m][1]) * (rk[1] - r[n][1])
                    + (rk[2] - r[m][2]) * (rk[2] - r[m][2]);

                    rkm = DIST(rk, r[m]);
                    rkm2 = SQUARE(rkm);

                    rkn = DIST(rk, r[n]);
                    rkn2 = SQUARE(rkn);

                    denom1 = rkm2 * (rkm - a);
                    denom2 = rkn2 * (rkn - a);
                    sum3 += rk_rm / (denom1 * denom2);
                }
            }
        }

        long double line4 = SQUARE(a) * sum3;

        return line1 + line2 + line3 + line4;
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
