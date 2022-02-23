#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include "rng.h"

#define ALPHA_INIT      0.1
#define ALPHA_UPDATE    0.1
#define N_ALPHA         10

#define MC_CYCLES       1000000L
#define EQUI_CYCLES     100000L
#define STEP_SIZE       0.5

#define N               100
#define D               3
#define OMEGA_HO        1
#define E_EXACT         D * 0.5 * N * OMEGA_HO

#define SAVE_RESULTS    false

double wave_function(double **r, double alpha) {
    double r2_sum = 0;
    for (int i = 0; i < N; ++i) {
        for (int dim = 0; dim < D; ++dim) {
            r2_sum += r[i][dim] * r[i][dim];
        }
    }
    return exp(- alpha * r2_sum);
}

double local_energy(double **r, double alpha) {
    double r2_sum = 0;
    for (int i = 0; i < N; ++i) {
        for (int dim = 0; dim < D; ++dim) {
            r2_sum += r[i][dim] * r[i][dim];
        }
    }
    return N * D * alpha - 2 * alpha * alpha * r2_sum + 0.5 * OMEGA_HO * OMEGA_HO * r2_sum;
}

void metropolis(double **output, unsigned int seed = ((unsigned) time(NULL))) {
    double **r_old = new double*[N];
    double **r_new = new double*[N];
    for (int i = 0; i < N; ++i) {
        r_old[i] = new double[D];
        r_new[i] = new double[D];
    }
    double wf_old, wf_new;
    double E_L;
    double E = 0, E2 = 0; 
    int idx_p, dim;

    printf("E_exact: %.3f \n", E_EXACT);
    RNG rng(seed);

    double alpha = ALPHA_INIT;
    for (int alpha_idx = 0; alpha_idx < N_ALPHA; ++alpha_idx) {
        auto start = std::chrono::steady_clock::now();

        for (idx_p = 0; idx_p < N; ++idx_p) {
            for (dim = 0; dim < D; ++dim) {
                r_old[idx_p][dim] = STEP_SIZE * (rng.rand_uniform() - 0.5);
            }   
        }
        wf_old = wave_function(r_old, alpha);

        for (long t = 0; t < MC_CYCLES; ++t) {
            for (idx_p = 0; idx_p < N; ++idx_p) {
                for (dim = 0; dim < D; ++dim) {
                    r_new[idx_p][dim] = r_old[idx_p][dim] + STEP_SIZE * (rng.rand_uniform() - 0.5);
                }
            }
            wf_new = wave_function(r_new, alpha);

            if (rng.rand_uniform() <= wf_new*wf_new / (wf_old*wf_old)) {
                for (idx_p = 0; idx_p < N; ++idx_p) {
                    for (dim = 0; dim < D; ++dim) {
                        r_old[idx_p][dim] = r_new[idx_p][dim];
                    }
                }
                wf_old = wf_new;
            }

            if (t > EQUI_CYCLES) {
                E_L = local_energy(r_old, alpha);
                E += E_L;
                E2 += E_L * E_L;
            }            
        }

        E /= (MC_CYCLES - EQUI_CYCLES);
        E2 /= (MC_CYCLES - EQUI_CYCLES);

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

    for (idx_p = 0; idx_p < N; ++idx_p) {
        delete[] r_new[idx_p];
        delete[] r_old[idx_p];
    }
    delete[] r_new;
    delete[] r_old;
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

