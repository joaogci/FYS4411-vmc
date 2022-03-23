#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include "../systems/system.h"
#include "../rng.h"

#include <cmath>

#define SQUARE(x)   ((x) * (x))

class MonteCarlo {
protected:

    System *system;
    RNG *rng;

    virtual long double acceptence_ratio(long double *r_new, int k) = 0;
    virtual void update_system(long double *r_new, int k) = 0;
    
public:

    MonteCarlo() {}

    ~MonteCarlo() {}

    long double step(long double *r_new, int k) {
        update_system(r_new, k);
        return acceptence_ratio(r_new, k);
    }

    void set_system(System *system_) {
        system = system_;
    }

    void set_rng(RNG *rng_) {
        rng = rng_;
    }

    virtual void print_info() = 0;
    virtual void set_param(double param_) = 0;

};

#endif // MONTE_CARLO_H
