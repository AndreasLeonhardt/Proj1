#include "trialfct.h"

// default constructor
TrialFct::TrialFct()
{
    funcParameters[0] = 1.0;
    funcParameters[1] = 1.0;
    stepwidth = 0.001;
    stepwidthsqr = stepwidth*stepwidth;
}

// constructor taking parameters from the config file
TrialFct::TrialFct(Config * parameters) : function(parameters)
{
    funcParameters[0] = parameters->lookup("a_min");
    funcParameters[1] = parameters->lookup("b_min");
    stepwidth = parameters->lookup("stepwidth");
    stepwidthsqr = stepwidth*stepwidth;
}


// calculate the value of the trial function the position given in the private variable position.
// Here the trial function is defined, currently
// f(r)=e^(-alpha(|r_1|+|r_2|) * e^(|r_1-r_2|/2*(1+beta|r_1-r_2|))
double TrialFct::getValue(positions * R)
{
    double result = 0.0;

    int nParticleshalf = nParticles/2;

    // create Slater matrix
    mat Slaterup  = zeros(nParticleshalf,nParticleshalf);
    mat Slaterdown  = zeros(nParticleshalf,nParticleshalf);



    for (int i=0;i<nParticleshalf;i++)
    {
        for (int j=0;j<nParticleshalf;j++)
        {
            Slaterup(i,j)   = hydrogen(j,               i,R);
            Slaterdown(i,j) = hydrogen(j+nParticleshalf,i,R);
        }
    }

    result = det(Slaterup)*det(Slaterdown);


    // jastrow factor

    double jastrow = 0.0;
    double dist;

    for (int i=0;i<nParticleshalf;i++)
    {   // equal spin, factor 1/2
        for (int j=0;j<i;j++)
        {
            dist = R->get_rr(i-1,j);
            jastrow += dist/(2+2*funcParameters[1]*dist);

            dist = R->get_rr(nParticleshalf+i-1,nParticleshalf+j);
            jastrow += dist/(2+2*funcParameters[1]*dist);
        }

        // opposite spin, factor 1/4
        for (int j=nParticleshalf;j<nParticles;j++)
        {
            dist = R->get_rr(j-1,i);
            jastrow += dist/(4+4*funcParameters[1]*dist);
        }

    }

    result *= exp(jastrow);

    return result;
}



// calculate the sum of the numerical second derivatives acting on the trail function ( \nabla^2_i f(x_1,...x_i,...x_n) ),
// derivatives act on the postion of particle n according to the given argument
// uses the position given in the privat variable position.
// using the numerical derivative: f''(x) = 1/h^2 * (f(x+h) + f(x-h) -2*f(x))
double TrialFct::getDivGrad(int particleNumber, positions * R)
{


    double value =-2*ndim*getValue(R);
    for (int i = 0; i<ndim; i++)
    {
        // move forward: r+h*e_i
        R->step(stepwidth,i,particleNumber);
        value += getValue(R);

        // move backwards: r-h*e_i
        R->step(-2*stepwidth,i,particleNumber);
        value += getValue(R);

        // move to middle
        R->step(stepwidth,i,particleNumber);

    }
    //cout << "h=" << stepwidth << "    h^2=" << stepwidthsqr << endl;
    return value/(stepwidthsqr);

}


// calculates numerically the quantum force, defined by 2/(f) * grad(f), where grad(f) is the gradient of f.
// particle number refers to the particle, on whichs position the derivatives act.
vec TrialFct::quantumForce(int particleNumber, positions *R)
{
    vec gradient = vec(ndim);
    for (int i =0; i<ndim; i++)
    {
        R->step(stepwidth,i,particleNumber);
        gradient(i)=getValue(R);

        R->step(-2*stepwidth,i,particleNumber);
        gradient(i)-=getValue(R);

        R->step(stepwidth,i,particleNumber);
    }

    gradient/=stepwidth*getValue(R);

    return gradient;

}


double TrialFct::hydrogen(int particleNumber, int orbital, positions * R)
{
    double result = 0;

    if (orbital==0)
    {
        result = exp(-funcParameters[0]*R->get_r(particleNumber));
    }

    else if(orbital==1)
    {
        result = ( 1.0 - 0.5* funcParameters[0] * R->get_r(particleNumber) )
                 * exp( -0.5* funcParameters[0] * R->get_r(particleNumber) );
    }
    else if(orbital==2 || orbital==3 || orbital==4)
    {
        result = funcParameters[0]*R->get_singlePos(particleNumber)(orbital-2)
                *exp(-0.5*funcParameters[0]*R->get_r(particleNumber));
    }

    return result;
}



// set and get private variables

//void TrialFct::setParameter(double newParameter, int parameterNumber)
//{
//    funcParameters[parameterNumber] = newParameter;
//}

//double TrialFct::getParameter(int parameterNumber)
//{
//    return funcParameters[parameterNumber];

//}


void TrialFct::set_stepwidth(double new_stepwidth)
{
    stepwidth = new_stepwidth;
    stepwidthsqr = stepwidth*stepwidth;

}

