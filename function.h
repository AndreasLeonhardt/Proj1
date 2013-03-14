#ifndef FUNCTION_H
#define FUNCTION_H


#include <armadillo>
#include <libconfig.h++>

#include "positions.h"

using namespace std;
using namespace arma;
using namespace libconfig;


class function
{
protected:

    // more general parameters (position vector, number of spacial dimensions, number of particles.
    int ndim, nParticles;
    double funcParameters[2];
    double stepwidth, stepwidthsqr;

public:
    function();
    function(Config * parameters);


    // calculate values of the function and its derivatives at the positon &R,
    // derivatives acting on the coordinates of particle particleNumber

    // value of the function
    virtual double getValue(positions * R)=0;

    // sum of the second derivatives, div(grad(f))
    virtual double getDivGrad(int particleNumber, positions * R)=0;

    // quantum force defined by g*grad(f)/f
    // g to be specified (1/2 for fermionic wave functions
    virtual vec quantumForce(int particleNumber, positions *R)=0;

    double get_stepwidth();
    double get_stepwidthsqr();

    void setParameter(double newParameter, int parameterNumber);
    double getParameter(int parameterNumber);


};

#endif // FUNCTION_H
