#include <iostream>

#include "system.h"
#include "rng.h"

System::System() {
    this->rng = new RNG((unsigned int) time(NULL));
}

System::System(unsigned int seed) {
    this->rng = new RNG(seed);
}

System::~System() {
    delete this->rng;

    if (this->r_allocated) {
        for (int i=0; i<N; i++) {
            delete[] this->r[i];
        }

        delete[] this->r;
        delete this->wf;
    }
}

void System::init_particles(int N, int dim, double alpha) {
    this->N = N;
    this->dim = dim;
    this->alpha = alpha;

    this->r = new long double*[N];
    for (int i=0; i<N; i++) {
        this->r[i] = new long double[dim];
    }

    this->wf = new WaveFunction(this->N, this->dim, this->alpha);

    this->r_allocated = true;
}

void System::set_simulation_params(long mc_cycles, double equi_fraction, double step_length) {
    this->mc_cycles = mc_cycles;
    this->equi_fraction = equi_fraction;
    this->step_length = step_length;

    this->set_params = true;
}



