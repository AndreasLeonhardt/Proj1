#ifndef MCINT_H
#define MCINT_H


#include <armadillo>
#include <libconfig.h++>

#include "trialfct.h"
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
    double stepSize;
    int acceptedSteps;
    double value;
    double value_analytical;
    double variance;
    double variance_analytical;

public:
    mcInt(Config * parameters);

    void integrate(TrialFct * fct, hamilton * H, Config *parameters);

    void set_nSamples(int NewnSamples);
    int  get_nSamples();

    void set_stepSize(double NewstepSize);
    int  get_stepSize();

    double get_acceptanceRatio();
    double get_value();
    double get_value_analytical();
    double get_variance();
    double get_variance_analytical();
};

#endif // MCINT_H
