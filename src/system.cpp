#include <iostream>
#include <cmath>
#include <chrono>
#include <vector>
#include <string>
#include <filesystem>
#include <fstream>
#include <iomanip>

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
    
    auto start = std::chrono::steady_clock::now();
    for (int t = 0; t < mc_cycles; t++) {
        idx_p = rng->rand() % N;
        wf_old = wave_function(r[idx_p]);

        for (d = 0; d < dim; d++) {
            r_new[d] = r[idx_p][d] + step_length * (rng->rand_uniform() - 0.5);
        }
        wf_new = wave_function(r_new);

        ratio = greens_function(r_new, r[idx_p]) * (wf_new*wf_new) / (wf_old*wf_old);
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
    
    auto end = std::chrono::steady_clock::now();
    double time =  (double) std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() * pow(10.0, -9);

    alphas.push_back(alpha);
    energies.push_back(E);
    energies2.push_back(E2);
    variances.push_back(E2 - E*E);
    wall_time_alpha.push_back(time);
    wall_time += time;

    printf("alpha: %.2f | E: %3.6f | var: %4.6f | time: %fs \n", 
                alphas.back(), 
                energies.back(), 
                variances.back(), 
                wall_time_alpha.back());
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

long double System::quantum_force(long double r) {
    return - 4 * alpha * r; 
}

long double System::greens_function(long double *r_new, long double *r_old) {
    long double gf = 0;

    for (int d = 0; d < dim; d++) {
        gf += 0.5 * (quantum_force(r_new[d]) - quantum_force(r_old[d])) * 
        (r_old[d] - r_new[d] + 0.5 * D * dt * (quantum_force(r_old[d] - quantum_force(r_new[d]))));
    }

    return exp(gf);
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

double System::get_wall_time() {
    return wall_time;
}

void System::write_results(std::string name, std::string path) {
    std::filesystem::create_directories(path); 

    std::ofstream file1(path + name);
    if (file1.is_open()) {
        file1 << "alpha,energy,energy2,variance,wall_time\n";
        for (int i = 0; i < alphas.size(); i++) {
            file1 << std::setprecision(10);
            file1 << alphas.at(i) << "," << energies.at(i) << "," << energies2.at(i) << "," << variances.at(i) << "," << wall_time_alpha.at(i) << "\n";
        }
        file1.close();
    } else {
        printf(" -- Error: can not open save file, please check you directory -- \n");
    }
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

