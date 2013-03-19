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


// returns 0 at the moment
double TrialFct_analytical::getDivGrad(int particleNumber, positions *R)
{
//    double value =-2*ndim*getValue(R);
//    for (int i = 0; i<ndim; i++)
//    {
//        // move forward: r+h*e_i
//        R->step(stepwidth,i,particleNumber);
//        value += getValue(R);

//        // move backwards: r-h*e_i
//        R->step(-2*stepwidth,i,particleNumber);
//        value += getValue(R);

//        // move to middle
//        R->step(stepwidth,i,particleNumber);
//    }
//    return value/(stepwidthsqr);
    cout << " analytical stuff"<<endl;
    return 0.0;
}

// returns the quantum force. well, 0 at the moment
vec TrialFct_analytical::quantumForce(int particleNumber, positions *R)
{
//    vec gradient = vec(ndim);
//    for (int i =0; i<ndim; i++)
//    {
//        R->step(stepwidth,i,particleNumber);
//        gradient(i)=getValue(R);

//        R->step(-2*stepwidth,i,particleNumber);
//        gradient(i)-=getValue(R);

//        R->step(stepwidth,i,particleNumber);
//    }

//    gradient/=stepwidth*getValue(R);
//    return gradient;
    return zeros(ndim);
}





