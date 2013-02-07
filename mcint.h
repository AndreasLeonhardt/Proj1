#ifndef MCINT_H
#define MCINT_H


#include <armadillo>
#include <libconfig.h++>

#include "trialfct.h"
#include "hamilton.h"

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
    double variance;

public:
    mcInt(Config * parameters);

    void integrate(TrialFct * fct, hamilton * H);

    void set_nSamples(int NewnSamples);
    int  get_nSamples();

    void set_stepSize(double NewstepSize);
    int  get_stepSize();

    double get_acceptanceRatio();
    double get_value();
    double get_variance();
};

#endif // MCINT_H
