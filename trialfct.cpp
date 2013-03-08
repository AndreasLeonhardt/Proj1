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
    double value = 0;

    for (int i=0;i<nParticles;i++)
    {
        value += -funcParameters[0]*R->get_r(i);
        for (int j=0;j<i;j++)
        {   double dist = R->get_rr(i-1,j);
            value += dist/(2+2*funcParameters[1]*dist);
        }
    }
    value = exp( value );

    return value;
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


// set and get private variables

void TrialFct::setParameter(double newParameter, int parameterNumber)
{
    funcParameters[parameterNumber] = newParameter;
}

double TrialFct::getParameter(int parameterNumber)
{
    return funcParameters[parameterNumber];

}


void TrialFct::set_stepwidth(double new_stepwidth)
{
    stepwidth = new_stepwidth;
    stepwidthsqr = stepwidth*stepwidth;

}

double TrialFct::get_stepwidth()
{
    return stepwidth;
}
