#include "trialfct_analytical.h"

TrialFct_analytical::TrialFct_analytical()
{
    alpha = 1.0;
    beta = 1.0;

}

TrialFct_analytical::TrialFct_analytical(Config * parameters)
{
    alpha = parameters->lookup("a_min");
    beta = parameters->lookup("b_min");
}


// calculate the value of the trial function the position given in the private variable position.
// Here the trial function is defined, currently
// f(r)=e^(-alpha(|r_1|+|r_2|) * e^(|r_1-r_2|/2*(1+beta|r_1-r_2|))
double TrialFct_analytical::getValue(positions * R)
{
    double value = 0;

    for (int i=0;i<nParticles;i++)
    {
        value += -alpha*R->get_r(i);
        for (int j=0;j<i;j++)
        {   double dist = R->get_rr(i-1,j);
            value += dist/(2+2*beta*dist);
        }
    }
    value = exp( value );

    return value;
}



double TrialFct_analytical::getDivGrad(int particleNumber, positions *R)
{
    return 0.0;
}

vec TrialFct_analytical::quantumForce(int particleNumber, positions *R)
{
    vec null =zeros(3);
    return null;
}


// set and get private variables

void TrialFct_analytical::set_alpha(double new_alpha)
{
    alpha = new_alpha;
}

double TrialFct_analytical::get_alpha()
{
    return alpha;
}


void TrialFct_analytical::set_beta(double new_beta)
{
    beta = new_beta;
}

double TrialFct_analytical::get_beta()
{
    return beta;
}
