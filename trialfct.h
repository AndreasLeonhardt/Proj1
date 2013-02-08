#ifndef TRIALFCT_H
#define TRIALFCT_H

#include <armadillo>
#include <libconfig.h++>

#include "positions.h"

using namespace std;
using namespace arma;
using namespace libconfig;


class TrialFct
{
    // parameters of the trial function
    double alpha, beta;

    // more general parameters (position vector, number of spacial dimensions, number of particles.
    int ndim, nParticles;



public:
    // constructor, takes the above parameters
    TrialFct(Config * parameters);

    // calculate the value of the trial function at a certain position
    double getValue(positions * R);


    // calculate the sum of the second derivatives acting on the trial function
    double getDivGrad(int particleNumber, positions * R);


    void set_alpha(double new_alpha);
    double get_alpha();

    void set_beta(double new_beta);
    double get_beta();


};

#endif // TRIALFCT_H

