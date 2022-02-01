#include <iostream>
#include <chrono>
#include <cmath>

#define ALPHA_INIT      0.2
#define ALPHA_UPDATE    0.05
#define N_ALPHA         10

#define MC_CYCLES       1000000L

#define N               100
#define OMEGA_HO        1
#define E_EXACT         0.5 * N * OMEGA_HO

#define RAND            (double) rand() / RAND_MAX

double wave_function(double *x, double alpha) {
    double x_sum2 = 0;
    for (int i = 0; i < N; ++i) {
        x_sum2 += x[i] * x[i];
    }
    return exp(- alpha * x_sum2);
}

double local_energy(double *x, double alpha) {
    double E_L = 0;
    for (int i = 0; i < N; ++i) {
        E_L += alpha - 2 * alpha * x[i] + 0.5 * OMEGA_HO * OMEGA_HO * x[i] * x[i];
    }
    return E_L;
}

void metropolis(double step_size, double **output, unsigned int seed = ((unsigned) time(NULL))) {
    double *x_old = new double[N];
    double *x_new = new double[N];
    double wf_old, wf_new;
    double delta_E;
    double E = 0, E2 = 0; 
    int idx_p;

    printf("E_exact: %.3f \n", E_EXACT);
    srand(seed);

    double alpha = ALPHA_INIT;
    for (int alpha_idx = 0; alpha_idx < N_ALPHA; ++alpha_idx) {
        auto start = std::chrono::steady_clock::now();

        for (idx_p = 0; idx_p < N; ++idx_p) {
            x_old[idx_p] = step_size * (RAND - 0.5);
        }
        wf_old = wave_function(x_old, alpha);

        for (long t = 0; t < MC_CYCLES; ++t) {
            for (idx_p = 0; idx_p < N; ++idx_p) {
                x_new[idx_p] = x_old[idx_p] + step_size * (RAND - 0.5);
            }
            wf_new = wave_function(x_new, alpha);

            if (RAND <= wf_new*wf_new / (wf_old*wf_old)) {
                for (idx_p = 0; idx_p < N; ++idx_p) {
                    x_old[idx_p] = x_new[idx_p];
                }
                wf_old = wf_new;
            }

            delta_E = local_energy(x_old, alpha);
            E += delta_E;
            E2 += delta_E * delta_E;
        }

        E /= MC_CYCLES;
        E2 /= MC_CYCLES;

        output[alpha_idx][0] = alpha;
        output[alpha_idx][1] = E;
        output[alpha_idx][2] = E2 - E*E;
        output[alpha_idx][3] = sqrt((E2 - E*E) / MC_CYCLES);
        
        auto end = std::chrono::steady_clock::now();
        double time = (double) std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() * pow(10.0, -9);
        printf("Alpha: %.2f | E: %.6f | Var: %.6f | Error: %.6f | Time: %.3fs \n", 
                output[alpha_idx][0], 
                output[alpha_idx][1], 
                output[alpha_idx][2], 
                output[alpha_idx][3],
                time);
        
        alpha += ALPHA_UPDATE;
    }

    delete[] x_new;
    delete[] x_old;
}

int main(int argc, char **argv) {
    double **output = new double *[N_ALPHA];
    for (int i = 0; i < N_ALPHA; ++i) {
        output[i] = new double[4];
    }

    auto start = std::chrono::steady_clock::now();
    metropolis(1.0, output);
    auto end = std::chrono::steady_clock::now();
    printf("Simulation time: %f seconds \n", (double) std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() * pow(10.0, -9));

    for (int i = 0; i < N_ALPHA; ++i) {
        delete[] output[i];
    }
    delete[] output;

    return 0;
}

