#include "positions.h"



positions::positions(Config * parameters)
{
    ndim = parameters->lookup("ndim");
    nParticles = parameters->lookup("nParticles");
    positions.randn(ndim,nParticles);
}



void positions::set_pos(mat * NewPositions)
{
    position = &NewPositions;

    for (int i=0;i<nParticles;i++)
    {
        // length of position vector
        r[i] = norm(position.col(i),2);

        // relative distance
        for (int j=0; j<i;j++)
        {
            rr(i-1,j) = norm(position.col(i)-position.col(j),2);
        }
    }
}

void positions::set_singlePos(vec * NewPosition, int particleNumber)
{
    // write new particle position
    position.col(particleNumber)=&NewPosition;

    // update r
    r[particleNumber]=norm(position.col(particleNumber),2);

    // update rr
    for (int j=0;j<particleNumber;j++)
    {
        rr(particleNumber-1,j)=norm(position.col(particleNumber)-position.col(j),2);
    }
    for (int i=particleNumber+1; i<nParticles;i++)
    {
        rr(i-1,particleNumber) = norm(position.col(particleNumber)-position.col(i),2);
    }
}


double positions::get_r(int i)
{
    return r[i];
}

double positions::get_rr(int i, int j)
{
    return rr(i,j);
}

mat positions::get_pos()
{
    return positions;
}


