#include <iostream>
#include <cmath>

#define ALPHA_INIT      0.4
#define ALPHA_UPDATE    0.05
#define N_ALPHA         20

#define MC_CYCLES       100000L

#define RAND            (double) rand() / RAND_MAX

double wave_function(double x, double alpha) {
    return exp(- 0.5 * alpha * alpha * x * x);
}

double local_energy(double x, double alpha) {
    return 0.5 * x * x * (1 - pow(alpha, 4)) + 0.5 * alpha * alpha;
}

void metropolis(double step_size, double **output, unsigned int seed = ((unsigned) time(NULL))) {
    double x_old, x_new;
    double wf_old, wf_new;
    double delta_E;
    double E = 0; 
    double E2 = 0;

    srand(seed);

    double alpha = ALPHA_INIT;
    for (int alpha_idx = 0; alpha_idx < N_ALPHA; ++alpha_idx) {
        x_old = step_size * (RAND - 0.5);
        wf_old = wave_function(x_old, alpha);

        for (long t = 0; t < MC_CYCLES; ++t) {
            x_new = x_old + step_size * (RAND - 0.5);
            wf_new = wave_function(x_new, alpha);

            if (RAND <= wf_new*wf_new / (wf_old*wf_old)) {
                x_old = x_new;
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

        printf("Alpha: %.2f | E: %.6f | Var: %.6f | Error: %.6f \n", 
                output[alpha_idx][0], 
                output[alpha_idx][1], 
                output[alpha_idx][2], 
                output[alpha_idx][3]);

        alpha += ALPHA_UPDATE;
    }
}

int main(int argc, char **argv) {
    double **output = new double *[N_ALPHA];
    for (int i = 0; i < N_ALPHA; ++i) {
        output[i] = new double[4];
    }

    metropolis(1.0, output);

    for (int i = 0; i < N_ALPHA; ++i) {
        delete[] output[i];
    }
    delete[] output;

    return 0;
}

