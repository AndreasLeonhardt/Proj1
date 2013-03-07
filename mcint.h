#ifndef MCINT_H
#define MCINT_H


#include <armadillo>
#include <libconfig.h++>

#include "function.h"
#include "hamilton.h"
#include "positions.h"
#include "lib.h"


using namespace std;
using namespace arma;
using namespace libconfig;


class mcInt
{
    int nParticles;
    int ndim;
    int nSamples;
    double timeStep;
    double sqrtTimeStep;
    int acceptedSteps;
    double value;
    double variance;

public:
    mcInt(Config * parameters);

    void integrate(function * fct, hamilton * H, Config *parameters, long idum);

    void set_nSamples(int NewnSamples);
    int  get_nSamples();

    void set_timeStep(double NewTimeStep);
    int  get_timeStep();

    double get_acceptanceRatio();
    double get_value();
    double get_variance();
};

#endif // MCINT_H
