#include "positions.h"



positions::positions(Config * parameters)
{
    ndim = parameters->lookup("ndim");
    nParticles = parameters->lookup("nParticles");
    pos.randn(ndim,nParticles);
    for (int i=0;i<nParticles;i++)
    {
        // length of position vector
        r[i] = norm(pos.col(i),2);

        // relative distance
        for (int j=0; j<i;j++)
        {
            rr(i-1,j) = norm(pos.col(i)-pos.col(j),2);
        }
    }
}



void positions::set_pos(mat  NewPositions)
{
    pos = NewPositions;

    for (int i=0;i<nParticles;i++)
    {
        // length of position vector
        r[i] = norm(pos.col(i),2);

        // relative distance
        for (int j=0; j<i;j++)
        {
            rr(i-1,j) = norm(pos.col(i)-pos.col(j),2);
        }
    }
}

void positions::set_singlePos(vec NewPosition, int particleNumber)
{
    // write new particle position
    pos.col(particleNumber)=NewPosition;

    // update r
    r[particleNumber]=norm(pos.col(particleNumber),2);

    // update rr
    for (int j=0;j<particleNumber;j++)
    {
        rr(particleNumber-1,j)=norm(pos.col(particleNumber)-pos.col(j),2);
    }
    for (int i=particleNumber+1; i<nParticles;i++)
    {
        rr(i-1,particleNumber) = norm(pos.col(particleNumber)-pos.col(i),2);
    }
}


void positions::step(double distance, int axis, int Particle)
{
    pos(axis,Particle)+=distance;
    // update r
    r[Particle]=norm(pos.col(Particle),2);

    // update rr
    for (int j=0;j<Particle;j++)
    {
        rr(Particle-1,j)=norm(pos.col(Particle)-pos.col(j),2);
    }
    for (int i=Particle+1; i<nParticles;i++)
    {
        rr(i-1,Particle) = norm(pos.col(Particle)-pos.col(i),2);
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
    return pos;
}


