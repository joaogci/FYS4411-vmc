#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include "../systems/system.h"
#include "../rng.h"

#include <cmath>

#define SQUARE(x)   x * x

class MonteCarlo {
protected:

    long mc_cycles;
    double equi_fraction;

    System *system;
    RNG *rng;
    
public:

    MonteCarlo(long mc_cycles_, double equi_fraction_, System *system_, RNG *rng_) {
        mc_cycles = mc_cycles_;
        equi_fraction = equi_fraction_;
        system = system_;
        rng = rng_;
    }

    ~MonteCarlo() {}

    long double acceptence_ratio(long double *r_new, long double *r_old);
    void update_system();

};

#endif // MONTE_CARLO_H
