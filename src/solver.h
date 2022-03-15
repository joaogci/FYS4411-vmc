#ifndef SOLVER_H
#define SOLVER_H

#include "rng.h"
#include "systems/system.h"
#include "samplers/monte_carlo.h"

class Solver {
private:

    RNG *rng;
    System* system;
    MonteCarlo* mc_sampler;

    long double *E;
    long double *E2;

public:

    Solver(System *system_, MonteCarlo* mc_sampler_, unsigned int seed_) {
        rng = new RNG(seed_);
        mc_sampler = mc_sampler_;
        system = system_;
    }
    
    ~Solver() {
        delete rng;
        delete system;
        delete mc_sampler;
    }

    void optimize_var_params();
    void solve();

    void write_results();
    

};

#endif // SOLVER_H
