#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include "../systems/system.h"
#include "../rng.h"

#include <cmath>

#define SQUARE(x)   (x * x)

class MonteCarlo {
protected:

    System *system;
    RNG *rng;

    virtual long double acceptence_ratio(long double *r_new, long double *r_old) = 0;
    virtual void update_system(long double *r_new, long double *r_old) = 0;
    
public:

    MonteCarlo() {}

    ~MonteCarlo() {}

    long double step(long double *r_new, long double *r_old) {
        update_system(r_new, r_old);
        return acceptence_ratio(r_new, r_old);
    }

    void set_system(System *system_) {
        system = system_;
    }

    void set_rng(RNG *rng_) {
        rng = rng_;
    }

};

#endif // MONTE_CARLO_H
