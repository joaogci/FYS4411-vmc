#ifndef SYSTEM_H
#define SYSTEM_H

#include <cmath>

#define SQUARE(x)       ((x) * (x))
#define DIST(rk, rm)    (sqrt(SQUARE(((rk[0]) - (rm[0]))) + SQUARE(((rk[1]) - (rm[1]))) + SQUARE(((rk[2]) - (rm[2])))))

class System {
protected:

    double a = 0.00433;         // Hard sphere radius
    double omega = 1.0; 
    double alpha = 0.0;         // Variational parameter
    double beta = sqrt(8.0);    // Another variational param that we don't vary

    bool check_init() {
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                if (DIST(r[i], r[j]) <= a) 
                    return true;
            }
        }
        return false;
    }

public:

    int N;                  // Number of particles
    int dim;                // Number of dimensions

    long double **r;        // Positions vector

    System(int N_, int dim_, double omega_) {
        N = N_;
        dim = dim_;
        omega = omega_;

        r = new long double*[N];
        do {
            for (int i = 0; i < N; i++) {
                r[i] = new long double[dim];

                for (int d = 0; d < dim; d++) {
                    r[i][d] = (((long double) rand() / RAND_MAX) - 0.5) * 2.0;
                }
            }
        } while (check_init());
    }

    System(int N_, int dim_, double omega_, double a_) {
        N = N_;
        dim = dim_;
        omega = omega_;
        a = a_;

        r = new long double*[N];
        do {
            for (int i = 0; i < N; i++) {
                r[i] = new long double[dim];

                for (int d = 0; d < dim; d++) {
                    r[i][d] = (((long double) rand() / RAND_MAX) - 0.5) * 2.0;
                }
            }
        } while (check_init());
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

    virtual long double evaluate_wf(long double *rk, int k) = 0;
    virtual void gradient_wf(long double *grad, long double *rk, int k) = 0;
    virtual long double laplacian_wf(long double *rk, int k) = 0;
    
    virtual void quantum_force(long double *force, long double *rk, int k) = 0;
    virtual long double local_energy() = 0;
    
    virtual void print_info() = 0;

};

#endif // SYSTEM_H
