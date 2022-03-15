#ifndef SYSTEM_H
#define SYSTEM_H

#include <cmath>

#define SQUARE(x)   x * x

class System {
protected:
    int N;
    int dim;

    // double a = 0.0;      // Hard sphere radius
    double omega;
    double alpha = 0;       // Variational parameter

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

    int get_dim() {
        return dim;
    }

    long double evaluate_wf(long double *r);
    void gradient_wf(long double *grad, long double *r);
    long double gradient_component_wf(long double x);
    long double laplacian_wf(long double *r);
    
    void quantum_force(long double *force, long double *r);
    long double local_energy();

};

#endif // SYSTEM_H
