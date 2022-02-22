#include <iostream>

#include "system.h"
#include "rng.h"

System::System() {
    rng = new RNG((unsigned int) time(NULL));
}

System::System(unsigned int seed) {
    rng = new RNG(seed);
}

void System::run_metropolis(double alpha_) {
    long double wf_old, wf_new;
    long double r_new[dim];
    long double ratio;
    int idx_p, d;

    long double E_l;
    long double E = 0, E2 = 0;

    alpha = alpha_;
    int counter = 0;
    
    for (int t = 0; t < mc_cycles; t++) {
        idx_p = rng->rand() % N;
        wf_old = wave_function(r[idx_p]);

        for (d = 0; d < dim; d++) {
            r_new[d] = r[idx_p][d] + step_length * (rng->rand_uniform() - 0.5);
        }
        wf_new = wave_function(r_new);

        ratio = (wf_new*wf_new) / (wf_old*wf_old);
        if (ratio >= 1 || rng->rand_uniform() <= ratio) {
            for (d = 0; d < dim; d++) {
                r[idx_p][d] = r_new[d];
            }
        }
        
        if (t >= mc_cycles * equi_fraction) {
            E_l = local_energy();
            E += E_l;
            E2 += E_l * E_l;
        }
    }

    E /= (mc_cycles * (1.0 - equi_fraction));
    E2 /= (mc_cycles * (1.0 - equi_fraction));

    printf("alpha: %.2f; E: %3.6f; var: %4.6f \n", alpha, E, E2 - E*E);
}

long double System::local_energy() {
    double r2_sum = 0;
    for (int i = 0; i < N; ++i) {
        for (int d = 0; d < dim; ++d) {
            r2_sum += r[i][d] * r[i][d];
        }
    }
    return N * dim * alpha - 2 * alpha * alpha * r2_sum + 0.5 * omega_ho * omega_ho * r2_sum;
}

long double System::wave_function() {
    double r2_sum = 0;
    for (int i = 0; i < N; ++i) {
        for (int d = 0; d < dim; ++d) {
            r2_sum += r[i][d] * r[i][d];
        }
    }
    return exp(- alpha * r2_sum);
}

long double System::wave_function(long double *r_i) {
    double r2_sum = 0;
    for (int d = 0; d < dim; ++d) {
        r2_sum += r_i[d] * r_i[d];
    }
    return exp(- alpha * r2_sum);
}

void System::init_particles(int N_, int dim_) {
    N = N_;
    dim = dim_;

    r = new long double*[N];
    for (int i = 0; i < N; i++) {
        r[i] = new long double[dim];
        for (int d = 0; d < dim; d++) {
            r[i][d] = equi_fraction * (rng->rand_uniform() - 0.5);
        }
    }

    r_allocated = true;
}

void System::init_hamiltonian(double omega_ho_) {
    omega_ho = omega_ho_;

    set_hamiltonian = true;
}

void System::set_simulation_params(long mc_cycles_, double equi_fraction_, double step_length_) {
    mc_cycles = mc_cycles_;
    equi_fraction = equi_fraction_;
    step_length = step_length_;

    set_params = true;
}

System::~System() {
    delete rng;

    if (r_allocated) {
        for (int i=0; i<N; i++) {
            delete[] r[i];
        }

        delete[] r;
    }
}

