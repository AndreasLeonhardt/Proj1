#include "function.h"

function::function()
{
    ndim = 3;
    nParticles = 2;
    nParticleshalf =1;
    nParams =2;
    R0=1.4;

    inverseSlaterDown = zeros(nParticles,nParticles);
    inverseSlaterUp   = zeros(nParticles,nParticles);

    for (int i=0;i<nParams;i++)
    {
        funcParameters[i]=1.0;
    }


}


function::function(Config * parameters)
{
    ndim = parameters->lookup("ndim");
    nParticles = parameters->lookup("nParticles");
    nParticleshalf = nParticles/2;
    nParams = parameters->lookup("nParameters");

    inverseSlaterDown = zeros(nParticles,nParticles);
    inverseSlaterUp   = zeros(nParticles,nParticles);

    funcParameters[0]=parameters->lookup("Parameters.[0]");
    funcParameters[1]=parameters->lookup("Parameters.[1]");
    R0=parameters->lookup("R0");
}





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

void function::setParameter(vec newParameters)
{
    for (int i =0;i<nParams;i++)
    {
        funcParameters[i] = newParameters[i];
    }
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


void function::set_R0(double R)
{
    R0=R;
}

double function::get_R0()
{
    return R0;
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
// modified in the lowest level to fit Hydrogen molecule.
double function::hydrogen(int particleNumber, int orbital, positions * R)
{
    double result;

    vec pos = R->get_singlePos(particleNumber);
    pos[2] -= R0/2;

    double r = norm(pos,2);

    result = exp(-funcParameters[0]*r);

    pos[2] +=R0;
    r=norm(pos,2);
    result += exp(-funcParameters[0]*r);


//    if (orbital==0)
//    {

//    }

//    else if(orbital==1)
//    {
//        result = ( 1.0 - 0.5* funcParameters[0] * r )
//                 * exp( -0.5* funcParameters[0] * r );
//    }
//    else if(orbital==2 || orbital==3 || orbital==4)
//    {
//        result = funcParameters[0]*R->get_singlePos(particleNumber)(orbital-2)
//                *exp(-0.5*funcParameters[0]*r);
//    }
//    else
//    {
//        cout << "single particle wave function undefined"<<endl;
//    }

    return result;
}

// gradient of hydrogen like wave functions
// hard coded as above
vec function::gradhydrogen(int particleNumber, int orbital, positions *R)
{
    vec result =zeros(ndim);
    vec pos = zeros(ndim);
    double a = funcParameters[0];

    pos = R->get_singlePos(particleNumber);
    pos[2]-=R0/2;
    double r = norm(pos,2);
    result =  -a/r *exp(-a*r)*pos;

    pos[2]+=R0;
    r = norm(pos,2);
    result +=  -a/r *exp(-a*r)*pos;

//    if (orbital==0)
//    {

//    }

//    else if(orbital==1)
//    {
//        result = R->get_singlePos(particleNumber);
//        result *= a/(4*r)*(a*r-4)*exp(-a*r/2);
//    }

//    else if(orbital==2 || orbital==3 || orbital==4)
//    {
//        result = R->get_singlePos(particleNumber)*a*R->get_singlePos(particleNumber)(orbital-2);
//        result(orbital-2)+= -2*r;
//        result *= -a/(2*r)*exp(-a*r/2);
//    }


    return result;
}


//laplace on hydrogen like wave functions
// closed form expressions, hard coded
double function::divgradhydrogen(int particleNumber, int orbital, positions* R)
{

    double result;
    double a = funcParameters[0];
    vec pos = R->get_singlePos(particleNumber);
    pos[2]-=R0/2;
    double r = norm(pos,2);
    result = a/r*(a*r-2) * exp(-a*r);

    pos[2]+=R0;
    r=norm(pos,2);
    result += a/r*(a*r-2) * exp(-a*r);


//    if (orbital==0)
//    {


//    }

//    else if(orbital==1)
//    {
//        result = -a/(8*r)*(a*a*r*r-10*a*r+18)*exp(-a*r/2);
//    }

//    else if(orbital==2 || orbital==3 || orbital==4)
//    {
//        vec RR = R->get_singlePos(particleNumber);
//        RR *= a*a/(4*r)*(a*r-8)*exp(-a*r/2);

//        result = RR(orbital-2);
//    }

    return result;

}



double function::dhydrogenda(int particleNumber, int orbital, positions *R)
{
    double result;
    double a = funcParameters[0];
    vec pos = R->get_singlePos(particleNumber);
    pos[2]-=R0/2;
    double r = norm(pos,2);
    result = -r * exp(-a*r);

    pos[2]+=R0;
    r=norm(pos,2);
    result += -r * exp(-a*r);

//    if (orbital==0)
//    {

//    }

//    else if(orbital==1)
//    {
//        result = (a/4*r*r-r)*exp(-a*r/2);
//    }


//    else if(orbital==2 || orbital==3 || orbital==4)
//    {
//        vec RR = R->get_singlePos(particleNumber);
//        RR *= (1-a*r/2)*exp(-a*r/2);

//        result = RR(orbital-2);
//    }

    return result;

}

