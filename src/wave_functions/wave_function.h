#ifndef WAVE_FUNCTION_H
#define WAVE_FUNCTION_H

#include <cmath>

class WaveFunction {
private:
    int N;
    int dim;
    double alpha;

public:
    WaveFunction(int N, int dim, double alpha) {
        this->N = N;
        this->dim = dim;
        this->alpha = alpha;
    }
    
    long double evaluate(long double **r) {
        double r2_sum = 0;
        for (int i=0; i<this->N; ++i) {
            for (int j=0; j<this->dim; ++j) {
                r2_sum += r[i][j] * r[i][j];
            }
        }

        return exp(- this->alpha * r2_sum);
    }

    long double laplacian(long double **r);

};

#endif
