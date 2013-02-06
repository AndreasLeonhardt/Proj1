#ifndef TRIALFCT_H
#define TRIALFCT_H

#include <armadillo>

using namespace std;
using namespace arma;


class TrialFct
{
    // parameters of the trial function
    double alpha, beta;

    // more general parameters (position vector, number of spacial dimensions, number of particles.
    mat position;
    int ndim, nParticles;


public:
    // constructor, takes the above parameters
    TrialFct(double a, double b, int d, int p);

    // calculate the value of the trial function at a certain position
    double getValue();

    // calculate the sum of the second derivatives acting on the trial function
    double getDivGrad(int particleNumber);

    // read and set the private values
    void set_position(mat newPosition);
    mat get_position();

    void set_alpha(double new_alpha);
    double get_alpha();

    void set_beta(double new_beta);
    double get_beta();

};

#endif // TRIALFCT_H

