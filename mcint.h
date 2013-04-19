#ifndef MCINT_H
#define MCINT_H


#include <armadillo>
#include <libconfig.h++>
#include <sys/types.h>
#include <sys/stat.h>


#include "function.h"
#include "hamilton.h"
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
    int SCnSamples;
    int SCthermalisationSteps;
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
    void integrate(function * fct, hamilton *H, long * idumadress, Config *parameters);
    vec StatGrad(function * fct, hamilton *H, long * idumadress, int nParams, Config *parameters);
    mat blocking(Config *parameters);

    void set_nSamples(int NewnSamples);
    int  get_nSamples();

    void set_timeStep(double NewTimeStep);
    int  get_timeStep();

    double get_acceptanceRatio();
    double get_value();
    double get_variance();
};

#endif // MCINT_H
