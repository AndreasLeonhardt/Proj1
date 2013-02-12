#ifndef HAMILTON_H
#define HAMILTON_H

#include <armadillo>
#include <libconfig.h++>

#include "positions.h"
#include "trialfct.h"

using namespace std;
using namespace arma;
using namespace libconfig;

class hamilton
{
    int ndim, nParticles, Z;

public:
    hamilton(Config * parameters);

    double localEnergy(TrialFct * trialfct, positions * R);


    double analytical_localEnergy(TrialFct * trialfct, positions * R);


};

#endif // HAMILTON_H