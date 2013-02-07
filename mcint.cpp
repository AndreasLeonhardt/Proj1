#include "mcint.h"

mcInt::mcInt(Config * parameters)
{
    nSamples = parameters->lookup("nSamples");
    stepSize = parameters->lookup("stepSize");
    ndim = parameters->lookup("ndim");
    nParticles = parameters->lookup("nParticles");
}




void mcInt::integrate(TrialFct * fct, hamilton * H)
{
    acceptedSteps = 0;
    value = 0.0;
    double stabw = 0.0;
    int i=0;

    double Phi = fct->getValue();
    P_old = Phi*Phi;
    double P_new;
    mat oldPosition(ndim,nParticles);
    mat newPosition(ndim,nParticles);

    while (i<nSamples)
    {

        // save old Position
        oldPosition = fct->get_position();
        R_old =fct->get_r();
        RR_old=fct->get_rr();

        // perform step
        newPosition = oldPosition+ (+2*randu<mat>(ndim,nParticles)-ones(ndim,nParticles))*stepSize;

        // calculate P_new/P_old
        fct->set_position(newPosition);
         Phi = fct->getValue();
         P_new =Phi*Phi;



        // compare to random variable
         vec zufallszahl(randu(1));
         if( zufallszahl(0) <= P_new/P_old )
        {


            // increase accepted steps
            acceptedSteps++;
        }
        else
        {
             fct->set_position(oldPosition);
        }

        // add energy value
         double ede = H->localEnergy(fct);
         value += ede;
         stabw += ede*ede;

        // increase i
        i++;

    }

    value /= nSamples*nParticles;
    variance = sqrt(stabw-value*value);
}




void mcInt::set_nSamples(int NewnSamples)
{
    nSamples = NewnSamples;
}

int  mcInt::get_nSamples()
{
    return nSamples;
}

void mcInt::set_stepSize(double NewstepSize)
{
    stepSize=NewstepSize;
}

int  mcInt::get_stepSize()
{
    return stepSize;
}


double mcInt::get_acceptanceRatio()
{
    return acceptedSteps/nSamples;
}

double mcInt::get_value()
{
    return value;
}


double mcInt::get_variance()
{
    return variance;
}
