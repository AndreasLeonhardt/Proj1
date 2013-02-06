#include "trialfct.h"

// constructor
TrialFct::TrialFct(double a, double b, int d, int p)
{
            alpha = a;
            beta = b;
            ndim = d;
            nParticles = p;


}


// calculate the value of the trial function the position given in the private variable position.
// Here the trial function is defined, currently
// f(r)=e^(-alpha(|r_1|+|r_2|) * e^(|r_1-r_2|/2*(1+beta|r_1-r_2|))
double TrialFct::getValue()
{
    double r[nParticles];
    for (int i=0; i<nParticles; i++)
    {
        r[i]=norm(position.col(i),2);
    }

    // This works just for 2 particles
    double r12= norm(position.col(0)-position.col(1),2);

    double value = exp( -alpha*(r[0]+r[1]) + r12/(2+2*beta*r12) );

    return value;
}

// calculate the sum of the numerical second derivatives acting on the trail function ( \nabla^2_i f(x_1,...x_i,...x_n) ),
// derivatives act on the postion of particle n according to the given argument
// uses the position given in the privat variable position.
double TrialFct::getDivGrad(int particleNumber)
{
    // stepwidth for numerical differenciation
    double h=0.001;

    double d=0.0;
    double step;

    for (int i = 0; i<ndim; i++)
    {   // r+h*e_i
        position(i,particleNumber)+=h;
        step = getValue();
        // r-h*e_i
        position(i,particleNumber)-=2*h;
        step += getValue();
        position(i,particleNumber)+=h;

        // f''(x)*h^2=f(x+h)+f(x-h) - 2*f(x)
        d+=step-2*getValue();
    }
    return d/(h*h);
}

// set and get private variables

void TrialFct::set_position(mat newPosition)
{
    position=newPosition;
}

mat TrialFct::get_position()
{
    return position;
}


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
