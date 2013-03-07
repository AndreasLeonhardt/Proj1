#include "hamilton.h"

// constructor
hamilton::hamilton(Config * parameters)
{
    ndim = parameters->lookup("ndim");
    nParticles = parameters->lookup("nParticles");
    Z = parameters->lookup("Z");
}

// calculate the local energy, 1/phi H phi.
// The Hamiltonian is defined within this function
// \hat{H} = \sum_n (\sum_i \frac12 \del_i^2) -\frac{Z}{r_n} +\sum_{j<i} \frac{1}{r_{ij}}

double hamilton::localEnergy(function * trialfct, positions * R)
{
}
