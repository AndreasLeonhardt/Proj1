#ifndef HAMILTON_ANALYTICAL_H
#define HAMILTON_ANALYTICAL_H


#include <armadillo>
#include <libconfig.h++>

#include "positions.h"
#include "function.h"
#include "hamilton.h"

using namespace std;
using namespace arma;
using namespace libconfig;


class hamilton_analytical : public hamilton
{
public:

    hamilton_analytical(Config * parameters);

    double localEnergy(function * trialfct, positions * R);

};

#endif // HAMILTON_ANALYTICAL_H