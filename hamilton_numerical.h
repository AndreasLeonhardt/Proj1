#ifndef HAMILTON_NUMERICAL_H
#define HAMILTON_NUMERICAL_H

#include <armadillo>
#include <libconfig.h++>

#include "positions.h"
#include "function.h"
#include "hamilton.h"

using namespace std;
using namespace arma;
using namespace libconfig;


class hamilton_numerical : public hamilton
{
public:
    hamilton_numerical(Config * parameters);
    double localEnergy(function *fct, positions * R);

};

#endif // HAMILTON_NUMERICAL_H
