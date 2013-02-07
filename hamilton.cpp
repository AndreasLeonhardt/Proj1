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

double hamilton::localEnergy(TrialFct * trialfct)
{


    double dens = 0;
    double value = trialfct->getValue();
            mat rr = trialfct->get_rr();

    for (int i=0;i<nParticles;i++)
    {
        // single particle energy
        // potential energy
        dens += -Z/trialfct->get_r()[i];
        // kinetic energy
        dens += -0.5*trialfct->getDivGrad(i)/value;

        // two paricle part
        for (int j=0; j<i;j++)
        {
            dens += 1/rr(i-1,j);
        }
    }

    return dens;
}
