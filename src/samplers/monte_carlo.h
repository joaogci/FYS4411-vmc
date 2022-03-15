#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include "../systems/system.h"
#include "../rng.h"

#include <cmath>

#define SQUARE(x)   x * x

class MonteCarlo {
protected:

    System *system;
    RNG *rng;
    
public:

    MonteCarlo(System *system_, RNG *rng_) {
        system = system_;
        rng = rng_;
    }

    ~MonteCarlo() {}

    long double acceptence_ratio(long double *r_new, long double *r_old);
    void update_system();

};

#endif // MONTE_CARLO_H
