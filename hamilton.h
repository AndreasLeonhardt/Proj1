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
protected:

    int ndim, nParticles;
    int Z[2];

public:

    hamilton(Config * parameters);

    double localEnergy(function *fct, positions * R);

};

#endif // HAMILTON_H
