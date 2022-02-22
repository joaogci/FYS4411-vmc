#ifndef SYSTEM_H
#define SYSTEM_H

#include <cmath>

#include "rng.h"

class System {
private:
    int N;                  // Number of particles
    int dim;                // Number of dimentions
    double alpha;           // Variational parameter
    
    long mc_cycles = 1000000L;         // Number of Monte Carlo (MC) cyles
    double equi_fraction = 0.3;   // Fraction of MC cycles used to achieve the equilibrium regime
    double step_length = 0.5;     // Step length for MC update
    
    double omega_ho = 1;        // Frequency of the HO 
    
    long double **r;        // (N * dim) array that contains the position of all particles

    RNG *rng;               // Random Number Generator

    bool r_allocated = false;
    bool set_params = false;
    bool set_hamiltonian = false;

public:
    System();
    System(unsigned int seed);
    ~System();

    void init_particles(int N_, int dim_);
    void set_simulation_params(long mc_cycles_, double equi_fraction_, double step_length_);
    void init_hamiltonian(double omega_ho_);
    long double local_energy();
    long double wave_function();
    long double wave_function(long double *r_i);
    void run_metropolis(double alpha_);
};

#endif // SYSTEM_H
