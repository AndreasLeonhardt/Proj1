#include "trialfct.h"

// constructor
TrialFct::TrialFct(Config * parameters)
{
    alpha = parameters->lookup("alpha");
    beta = parameters->lookup("beta");
            ndim = parameters->lookup("ndim");
            nParticles = parameters->lookup("nParticles");

}


// calculate the value of the trial function the position given in the private variable position.
// Here the trial function is defined, currently
// f(r)=e^(-alpha(|r_1|+|r_2|) * e^(|r_1-r_2|/2*(1+beta|r_1-r_2|))
double TrialFct::getValue(positions * R)
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



// calculate the sum of the numerical second derivatives acting on the trail function ( \nabla^2_i f(x_1,...x_i,...x_n) ),
// derivatives act on the postion of particle n according to the given argument
// uses the position given in the privat variable position.
double TrialFct::getDivGrad(int particleNumber, positions * R)
{
    // stepwidth for numerical differenciation
    double h=0.01;

    double value;

    for (int i = 0; i<ndim; i++)
    {
        // move forward: r+h*e_i

        R->step(h,i,particleNumber);
        value = getValue(R);


        // move backwards: r-h*e_i
        R->step(-2*h,i,particleNumber);
        value += getValue(R);

        // move to middle
        R->step(h,i,particleNumber);
        value+=-2*getValue(R);
    }
    return value/(h*h);
}



// set and get private variables

void TrialFct::set_alpha(double new_alpha)
{
    alpha = new_alpha;
}

double TrialFct::get_alpha()
{
    return alpha;
}


void TrialFct::set_beta(double new_beta)
{
    beta = new_beta;
}

double TrialFct::get_beta()
{
    return beta;
}
