#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <string>
#include <filesystem>
#include <fstream>
#include <chrono>

#include "rng.h"
#include "systems/system.h"
#include "samplers/monte_carlo.h"

class Solver {
private:

    RNG *rng;
    System* system;
    MonteCarlo* mc_sampler;

    double alpha;
    long double *E_sampled = NULL;
    long double *E2_sampled = NULL;
    long double run_time = 0;
    double accepted_ratio;

    long mc_cycles;
    double equi_fraction;
    double measure_after;
    long measure_cycles;
    int print_skip;

    double eta = 1e-4;
    double tol = 1e-6;
    double h = 1e-3;

    double solve_optimizer(double alpha_) {
        int idx_p;
        long double r_new[system->dim];
        long double ratio;

        long double E = 0;

        system->set_var_params(alpha_);

        for (int t = 0; t < mc_cycles; t++) {
            idx_p = rng->rand() % system->N;

            ratio = mc_sampler->step(r_new, idx_p);
            if (ratio >= 1 || rng->rand_uniform() <= ratio) {
                for (int d = 0; d < system->dim; d++) {
                    system->r[idx_p][d] = r_new[d];
                }
            }

            if(t >= measure_after) {
                E += system->local_energy();
            }
        }

        return E / measure_cycles;
    }

public:

    Solver(System *system_, MonteCarlo* mc_sampler_, long mc_cycles_, double equi_fraction_, unsigned int seed_) {
        rng = new RNG(seed_);
        mc_sampler = mc_sampler_;
        system = system_;

        mc_sampler->set_rng(rng);
        mc_sampler->set_system(system);

        mc_cycles = mc_cycles_;
        equi_fraction = equi_fraction_;
        measure_after = mc_cycles * equi_fraction;
        measure_cycles = mc_cycles * (1.0 - equi_fraction);

        E_sampled = new long double[measure_cycles];
        E2_sampled = new long double[measure_cycles];
    }

    ~Solver() {
        delete rng;
        delete system;
        delete mc_sampler;
    }

    void set_solve_params(long mc_cycles_, double equi_fraction_) {
        mc_cycles = mc_cycles_;
        equi_fraction = equi_fraction_;
        measure_after = mc_cycles * equi_fraction;
        measure_cycles = mc_cycles * (1.0 - equi_fraction);
        print_skip = mc_cycles / 2;

        if (E_sampled != NULL) {
            delete[] E_sampled;
        }
        if (E2_sampled != NULL) {
            delete[] E2_sampled;
        }

        E_sampled = new long double[mc_cycles];
        E2_sampled = new long double[mc_cycles];
    }

    void solve(double alpha_) {
        int idx_p;
        long double r_new[system->dim];
        long double ratio;

        long double E_l, E = 0, E2 = 0;

        alpha = alpha_;
        system->set_var_params(alpha);

        printf("Starting VMC calculation with \nalpha: %.4f | mc_cycles: 2^%i | measure_after: 2^%i \n", alpha, (int) log2(mc_cycles), (int) log2(measure_after));
        system->print_info();
        mc_sampler->print_info();
        printf("\n");

        auto start = std::chrono::steady_clock::now();
        for (long t = 0; t < mc_cycles; t++) {
            idx_p = rng->rand() % system->N;

            ratio = mc_sampler->step(r_new, idx_p);
            if (ratio >= 1 || rng->rand_uniform() <= ratio) {
                for (int d = 0; d < system->dim; d++) {
                    system->r[idx_p][d] = r_new[d];
                }
                accepted_ratio++;
            }

            E_l = system->local_energy();
            E_sampled[t] = E_l;
            E2_sampled[t] = E_l * E_l;

            if(t >= measure_after) {
                E += E_sampled[t];
                E2 += E2_sampled[t];
            }

            if (t == print_skip) {
                auto end = std::chrono::steady_clock::now();
                run_time =  (double) std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() * pow(10.0, -9);
                printf("Half way through the computation in %Lfs \n", run_time);
            }
        }
        auto end = std::chrono::steady_clock::now();
        run_time =  (double) std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() * pow(10.0, -9);

        E /= measure_cycles;
        E2 /= measure_cycles;
        accepted_ratio /= mc_cycles;

        printf("alpha: %.4f | E: %4.6Lf | var: %4.6Lf | time: %Lfs \n\n", alpha_, E, E2 - E*E, run_time);
    }

    void set_optimization_params(double eta_, double tol_, double h_, long mc_cycles_, double equi_fraction_) {
        mc_cycles = mc_cycles_;
        equi_fraction = equi_fraction_;
        measure_after = mc_cycles * equi_fraction;
        measure_cycles = mc_cycles * (1.0 - equi_fraction); 

        eta = eta_;
        tol = tol_;
        h = h_;
    } 

    double optimize_var_params(double alpha_0_) {
        double opt_alpha = alpha_0_;
        long double derivative = tol + 1;
        long double Ep, Em;

        long count = 0;

        printf("Starting optimization with \neta: %.7f | tol: %.6f | h: %.6f | mc_cycles 2^%i \n", eta, tol, h, (int) log2(mc_cycles));
        system->print_info();
        printf("\n");

        auto start = std::chrono::steady_clock::now();
        while (fabs(derivative) >= tol) {
            Ep = solve_optimizer(opt_alpha + h);
            Em = solve_optimizer(opt_alpha - h);

            derivative = eta * (Ep - Em) / (2.0 * h);
            opt_alpha += - derivative;

            if (count % 10 == 0) {
                printf("iter: %5li | alpha: %.4f | delta_alpha: %6.5Lf \n", count, opt_alpha, derivative);
            }
            count++;
        }
        auto end = std::chrono::steady_clock::now();
        double time =  (double) std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() * pow(10.0, -9);

        printf("opt_alpha: %.4f  | iter: %5li | time: %f\n\n", opt_alpha, count, time);

        return opt_alpha;
    }

    void write_results(std::string name, std::string path) {
        std::string file_loc = path + name;
        std::filesystem::create_directories(path);

        std::ofstream file1(file_loc);
        if (file1.is_open()) {
            file1 << std::setprecision(10);
            file1 << "mc_cycles,measure_cycles,alpha,runtime,acceptance_ratio\n";
            file1 << mc_cycles << "," << measure_cycles << "," << alpha << "," << run_time << "," << accepted_ratio << "\n";
            file1 << "energy,energy2\n";
            for (int i = 0; i < mc_cycles; i++) {
                file1 << E_sampled[i] << "," << E2_sampled[i] << "\n";
            }
            file1.close();
        } else {
            printf(" -- Error: can not open save file: \"%s\" -- \n", file_loc.c_str());
        }

        printf("Written results to file \"%s\" with success \n", file_loc.c_str());
    }

};

#endif // SOLVER_H
