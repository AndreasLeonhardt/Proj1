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

    int ndim, nParticles, Z;

public:

    hamilton(Config * parameters);

    virtual double localEnergy(function *fct, positions * R)=0;

};

#endif // HAMILTON_H
