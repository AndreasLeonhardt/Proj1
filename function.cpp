#include "function.h"

function::function()
{
    ndim = 3;
    nParticles = 2;
}


function::function(Config * parameters)
{
    ndim = parameters->lookup("ndim");
    nParticles = parameters->lookup("nParticles");
}


void function::setParameter(double newParameter, int parameterNumber)
{
}

double function::getParameter(int parameterNumber)
{
}


// calculate values of the function and its derivatives at the positon &R,
// derivatives acting on the coordinates of particle particleNumber

// value of the function
double function::getValue(positions * R)
{}

// sum of the second derivatives, div(grad(f))
double function::getDivGrad(int particleNumber, positions * R)
{}

// quantum force defined by g*grad(f)/f
// g to be specified (1/2 for fermionic wave functions
vec function::quantumForce(int particleNumber, positions *R)
{
}
