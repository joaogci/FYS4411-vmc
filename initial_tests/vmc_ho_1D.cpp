#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>

#define ALPHA_INIT      0.05
#define ALPHA_UPDATE    0.05
#define N_ALPHA         20

#define MC_CYCLES       1000000L
#define STEP_SIZE       1.0

#define N               100
#define OMEGA_HO        1
#define E_EXACT         0.5 * N * OMEGA_HO

#define RAND            (double) rand() / RAND_MAX
#define SAVE_RESULTS    false

double wave_function(double *x, double alpha) {
    double x2_sum = 0;
    for (int i = 0; i < N; ++i) {
        x2_sum += x[i] * x[i];
    }
    return exp(- alpha * x2_sum);
}

double local_energy(double *x, double alpha) {
    double E_L = 0;
     for (int i = 0; i < N; ++i) {
        E_L += alpha - 2 * alpha * alpha * x[i] * x[i] + 0.5 
                * OMEGA_HO * OMEGA_HO * x[i] * x[i];
    }
    return E_L;
}

void metropolis(double **output, unsigned int seed = ((unsigned) time(NULL))) {
    double *x_old = new double[N];
    double *x_new = new double[N];
    double wf_old, wf_new;
    double E_L;
    double E = 0, E2 = 0; 
    int idx_p;

    printf("E_exact: %.3f \n", E_EXACT);
    srand(seed);

    double alpha = ALPHA_INIT;
    for (int alpha_idx = 0; alpha_idx < N_ALPHA; ++alpha_idx) {
        auto start = std::chrono::steady_clock::now();

        for (idx_p = 0; idx_p < N; ++idx_p) {
            x_old[idx_p] = STEP_SIZE * (RAND - 0.5);
        }
        wf_old = wave_function(x_old, alpha);

        for (long t = 0; t < MC_CYCLES; ++t) {
            for (idx_p = 0; idx_p < N; ++idx_p) {
                x_new[idx_p] = x_old[idx_p] + STEP_SIZE * (RAND - 0.5);
            }
            wf_new = wave_function(x_new, alpha);

            if (RAND <= wf_new*wf_new / (wf_old*wf_old)) {
                for (idx_p = 0; idx_p < N; ++idx_p) {
                    x_old[idx_p] = x_new[idx_p];
                }
                wf_old = wf_new;
            }

            E_L = local_energy(x_old, alpha);
            E += E_L;
            E2 += E_L * E_L;
        }

        E /= MC_CYCLES;
        E2 /= MC_CYCLES;

        output[alpha_idx][0] = alpha;
        output[alpha_idx][1] = E;
        output[alpha_idx][2] = E2 - E*E;
        
        auto end = std::chrono::steady_clock::now();
        double time = (double) std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() * pow(10.0, -9);
        printf("Alpha: %.2f | E: %.6f | E2: %.6f | Var: %.6f | Time: %.3fs \n", 
                output[alpha_idx][0], 
                output[alpha_idx][1], 
                E2, 
                output[alpha_idx][2],
                time);
        
        alpha += ALPHA_UPDATE;
    }

    delete[] x_new;
    delete[] x_old;
}

int main(int argc, char **argv) {
    double **output = new double *[N_ALPHA];
    for (int i = 0; i < N_ALPHA; ++i) {
        output[i] = new double[3];
    }

    auto start = std::chrono::steady_clock::now();
    metropolis(output);
    auto end = std::chrono::steady_clock::now();
    printf("Simulation time: %f seconds \n", (double) std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() * pow(10.0, -9));

    if (SAVE_RESULTS) {
        std::ofstream file("results/vmc_ho_N_" + std::to_string(N) + "_omega_" + std::to_string(OMEGA_HO) + "_MC_1E" + std::to_string((int) log10(MC_CYCLES)) + "_na_" + std::to_string(N_ALPHA) + ".txt");
        if (file.is_open()) {
            for (int i = 0; i < N_ALPHA; ++i) {
                file << output[i][0] << "," << output[i][1] << "," << output[i][2] << "\n";
            }
            file.close();
        }
    }
    

    for (int i = 0; i < N_ALPHA; ++i) {
        delete[] output[i];
    }
    delete[] output;

    return 0;
}

