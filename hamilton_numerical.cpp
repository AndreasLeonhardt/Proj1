#include "hamilton_numerical.h"



hamilton_numerical::hamilton_numerical(Config * parameters):hamilton(parameters)
{
}


double hamilton_numerical::localEnergy(function *fct, positions * R)
{

    double dens = 0;
    double value = fct->getValue(R);

    for (int i=0;i<nParticles;i++)
    {
        // single particle energy
        // potential energy
        dens += -Z/(R->get_r(i));



        // kinetic energy
        dens += -0.5*(fct->getDivGrad(i,R))/value;


        // two paricle part
        for (int j=0; j<i;j++)
        {
            dens += 1/R->get_rr(i-1,j);

        }
    }

    return dens;
}

