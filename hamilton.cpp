#include "hamilton.h"

// constructor
hamilton::hamilton(Config * parameters)
{
    ndim = parameters->lookup("ndim");
    nParticles = parameters->lookup("nParticles");
    Z = parameters->lookup("Z");
}


double hamilton::localEnergy(function *fct, positions *R)
{
    double dens = 0.0;

    for (int i=0;i<nParticles;i++)
    {
        // single particle energy
        // potential energy
        dens += -Z/(R->get_r(i));



        // kinetic energy
        dens += -0.5*(fct->getDivGradOverFct(i,R));
     //   cout <<"kinetic energy = "<<fct->getDivGradOverFct(i,R)<<endl;

        // two paricle part
        for (int j=0; j<i;j++)
        {
            dens += 1/R->get_rr(i-1,j);

        }
    }

    return dens;
}

