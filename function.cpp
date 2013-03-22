#include "function.h"

function::function()
{
    ndim = 3;
    nParticles = 2;
    nParticleshalf =1;

    inverseSlaterDown = zeros(nParticles,nParticles);
    inverseSlaterUp   = zeros(nParticles,nParticles);


}


function::function(Config * parameters)
{
    ndim = parameters->lookup("ndim");
    nParticles = parameters->lookup("nParticles");
    nParticleshalf = nParticles/2;

    inverseSlaterDown = zeros(nParticles,nParticles);
    inverseSlaterUp   = zeros(nParticles,nParticles);

}




// calculate values of the function and its derivatives at the positon &R,
// derivatives acting on the coordinates of particle particleNumber

//// value of the function
//double function::getValue(positions * R)
//{}

//// sum of the second derivatives, div(grad(f))
//double function::getDivGrad(int particleNumber, positions * R)
//{}

//// quantum force defined by g*grad(f)/f
//// g to be specified (1/2 for fermionic wave functions
//vec function::quantumForce(int particleNumber, positions *R)
//{
//}


double function::get_stepwidth()
{
    return stepwidth;
}

double function::get_stepwidthsqr()
{
    return stepwidth;
}


void function::setParameter(double newParameter, int parameterNumber)
{
    funcParameters[parameterNumber]=newParameter;
}

double function::getParameter(int parameterNumber)
{
    return funcParameters[parameterNumber];
}


void function::setndim(int newvalue)
{
    ndim =newvalue;
}

int function::getndim()
{
    return ndim;
}

void function::setnParticles(int numberofParticles)
{
    nParticles=numberofParticles;
    nParticleshalf =nParticles/2;
}

int function::getnParticles()
{
    return nParticles;
}

mat function::getinvslatermatrix(int i)
{   if(i)
    {
        return inverseSlaterUp;
    }
    else
    {
        return inverseSlaterDown;
    }
}
