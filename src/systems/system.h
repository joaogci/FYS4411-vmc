#ifndef SYSTEM_H
#define SYSTEM_H

#include <cmath>

#define SQUARE(x)       (x * x)
#define DIST(rk, rm)    (sqrt(SQUARE(rk[0] - rm[0]) + SQUARE(rk[1] - rm[1]) + SQUARE(rk[2] - rm[2])))

class System {
protected:
    int N;
    int dim;

    double a = 0.0;         // Hard sphere radius
    double omega;
    double alpha = 0.0;     // Variational parameter
    double beta = 1.0;      // Another variational param that we don't vary

    long double **r;        // Positions vector

public:

    System(int N_, int dim_, double omega_) {
        N = N_;
        dim = dim_;
        omega = omega_;

        r = new long double*[N];
        for (int i = 0; i < N; i++) {
            r[i] = new long double[dim];

            for (int d = 0; d < dim; d++) {
                r[i][d] = ((long double) rand() / RAND_MAX) - 0.5;
            }
        }


    }
    ~System() {
        for (int i = 0; i < N; i++) {
            delete[] r[i];
        }
        delete[] r;
    }

    void set_var_params(double alpha_) {
        alpha = alpha_;
    }

    int get_N() {
        return N;
    }

    int get_dim() {
        return dim;
    }

    long double *get_rk(int k) {
        return r[k];
    }

    void update_rk(int k, long double *r_new) {
        for (int d = 0; d < dim; d++) {
            r[k][d] = r_new[d];
        }
    }

    virtual long double evaluate_wf(long double *r) = 0;
    virtual void gradient_wf(long double *grad, long double *r) = 0;
    virtual long double gradient_component_wf(long double x) = 0;
    virtual long double laplacian_wf(long double *rk, int k) = 0;
    
    virtual void quantum_force(long double *force, long double *r) = 0;
    virtual long double local_energy() = 0;

};

#endif // SYSTEM_H
