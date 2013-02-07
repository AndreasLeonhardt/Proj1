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
    mat position;
    int ndim, nParticles;
    vec r;
    mat rr;


public:
    // constructor, takes the above parameters
    TrialFct(Config * parameters, double a, double b );

    // calculate the value of the trial function at a certain position
    double getValue();
    double getValue(mat R,mat RR);

    // calculate the sum of the second derivatives acting on the trial function
    double getDivGrad(int particleNumber);

    // read and set the private values
    void set_position(mat newPosition);
    mat get_position();

    void updateParticelPosition(vec newPosition,int particleNumber);

    void set_alpha(double new_alpha);
    double get_alpha();

    void set_beta(double new_beta);
    double get_beta();

    vec get_r();
    mat get_rr();

};

#endif // TRIALFCT_H

