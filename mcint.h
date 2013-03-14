#ifndef MCINT_H
#define MCINT_H


#include <armadillo>
#include <libconfig.h++>

#include "function.h"
#include "hamilton.h"
#include "hamilton_numerical.h"
#include "positions.h"



using namespace std;
using namespace arma;
using namespace libconfig;


class mcInt
{
    int nParticles;
    int ndim;
    int nSamples;
    int thermalisationSteps;
    double timeStep;
    double sqrtTimeStep;
    int acceptedSteps;
    double value;
    double variance;
    double P_old;


public:
    mcInt();
    mcInt(Config * parameters);

    positions * Step(function * fct, positions *Rold, long * idumadress, Config *parameters);
    positions * thermalise(function * fct, long * idumadress, Config *parameters);
    void integrate(function * fct, hamilton *H, positions * Rold, long * idumadress, Config *parameters);

    void set_nSamples(int NewnSamples);
    int  get_nSamples();

    void set_timeStep(double NewTimeStep);
    int  get_timeStep();

    double get_acceptanceRatio();
    double get_value();
    double get_variance();
};

#endif // MCINT_H
