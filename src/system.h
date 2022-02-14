#ifndef SYSTEM_H
#define SYSTEM_H

#include "rng.h"
#include "sampler.h"
#include "wave_functions/wave_function.h"

class System {
private:
    int N;                  // Number of particles
    int dim;                // Number of dimentions
    long mc_cycles;         // Number of Monte Carlo (MC) cyles
    double equi_fraction;   // Fraction of MC cycles used to achieve the equilibrium regime
    double step_length;     // Step length for MC update
    double alpha;           // Variational parameter
    
    long double **r;        // N * dim array that contains the position of all particles

    RNG *rng;               // Random Number Generator
    Sampler *sampler;       // Class Sampler instance
    WaveFunction *wf;       // Trial WaveFunction for the system
    

    bool r_allocated = false;
    bool set_params = false;

public:
    System();
    System(unsigned int seed);
    ~System();

    void init_particles(int N, int dim, double alpha);
    void set_simulation_params(long mc_cycles, double equi_fraction, double step_length);

};

#endif // SYSTEM_H
