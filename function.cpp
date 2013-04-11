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

mat function::getinvslatermatrix(int particleNumber)
{   if(particleNumber<nParticleshalf)
    {
        return inverseSlaterUp;
    }
    else
    {
        return inverseSlaterDown;
    }
}

void function::setSlaterinv(int particleNumber, mat newSlaterInv)
{
    if (particleNumber<nParticleshalf)
    {
        inverseSlaterUp = newSlaterInv;
    }
    else
    {
        inverseSlaterDown = newSlaterInv;
    }
}

//===========================================================================================================
// hydrogen like wave functions, with parameter alpha
// hard coded for the first levels up to n=2, l=1, m_l = -1,0,1
// not including spin degeneracy.
double function::hydrogen(int particleNumber, int orbital, positions * R)
{
    double r=R->get_r(particleNumber);

    double result;

    if (orbital==0)
    {
        result = exp(-funcParameters[0]*r);
    }

    else if(orbital==1)
    {
        result = ( 1.0 - 0.5* funcParameters[0] * r )
                 * exp( -0.5* funcParameters[0] * r );
    }
    else if(orbital==2 || orbital==3 || orbital==4)
    {
        result = funcParameters[0]*R->get_singlePos(particleNumber)(orbital-2)
                *exp(-0.5*funcParameters[0]*r);
    }
    else
    {
        cout << "single particle wave function undefined"<<endl;
    }

    return result;
}

// gradient of hydrogen like wave functions
// hard coded as above
vec function::gradhydrogen(int particleNumber, int orbital, positions *R)
{
    vec result =zeros(ndim);
    double a = funcParameters[0];
    double r = R->get_r(particleNumber);

    if (orbital==0)
    {
        result = R->get_singlePos(particleNumber);
        result *=  -a/r *exp(-a*r);
    }

    else if(orbital==1)
    {
        result = R->get_singlePos(particleNumber);
        result *= a/(4*r)*(a*r-4)*exp(-a*r/2);
    }

    else if(orbital==2 || orbital==3 || orbital==4)
    {
        result = R->get_singlePos(particleNumber)*a*R->get_singlePos(particleNumber)(orbital-2);
        result(orbital-2)+= -2*r;
        result *= -a/(2*r)*exp(-a*r/2);
    }


    return result;
}


//laplace on hydrogen like wave functions
// closed form expressions, hard coded
double function::divgradhydrogen(int particleNumber, int orbital, positions* R)
{

    double result;
    double a = funcParameters[0];
    double r = R->get_r(particleNumber);

    if (orbital==0)
    {

        result = a/r*(a*r-2) * exp(-a*r);
    }

    else if(orbital==1)
    {
        result = -a/(8*r)*(a*a*r*r-10*a*r+18)*exp(-a*r/2);
    }

    else if(orbital==2 || orbital==3 || orbital==4)
    {
        vec RR = R->get_singlePos(particleNumber);
        RR *= a*a/(4*r)*(a*r-8)*exp(-a*r/2);

        result = RR(orbital-2);
    }

    return result;

}



double function::dhydrogenda(int particleNumber, int orbital, positions *R)
{
    double result;
    double a = funcParameters[0];
    double r = R->get_r(particleNumber);

    if (orbital==0)
    {
        result = -r * exp(-a*r);
    }

    else if(orbital==1)
    {
        result = (a/4*r*r-r)*exp(-a*r/2);
    }


    else if(orbital==2 || orbital==3 || orbital==4)
    {
        vec RR = R->get_singlePos(particleNumber);
        RR *= (1-a*r/2)*exp(-a*r/2);

        result = RR(orbital-2);
    }

    return result;

}

