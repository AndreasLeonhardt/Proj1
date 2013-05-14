#include "hamilton.h"

// constructor
hamilton::hamilton(Config * parameters)
{
    ndim = parameters->lookup("ndim");
    nParticles = parameters->lookup("nParticles");
    Z[0] = parameters->lookup("Z.[0]");
    Z[1] = parameters->lookup("Z.[1]");
}


double hamilton::localEnergy(function *fct, positions *R)
{
    double dens = 0.0;
    double r;
    double R0=fct->get_R0();
    dens +=1/R0;

    for (int i=0;i<nParticles;i++)
    {
        // single particle energy
        // potential energy
        vec pos = R->get_singlePos(i);
        pos[2]-=R0/2;
        r=norm(pos,2);
        dens += -Z[0]/r;

        pos[2]+=R0;
        r=norm(pos,2);
        dens += -Z[1]/r;


        // kinetic energy
        dens += -0.5*(fct->getDivGradOverFct(i,R));

        // two paricle part
        for (int j=0; j<i;j++)
        {
            dens += 1/R->get_rr(i-1,j);
        }
    }

    return dens;
}

