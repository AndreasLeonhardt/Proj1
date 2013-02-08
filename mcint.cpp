#include "mcint.h"

mcInt::mcInt(Config * parameters)
{
    nSamples = parameters->lookup("nSamples");
    stepSize = parameters->lookup("stepSize");
    ndim = parameters->lookup("ndim");
    nParticles = parameters->lookup("nParticles");
}




void mcInt::integrate(TrialFct * fct, hamilton * H, Config * parameters)
{
    acceptedSteps = 0;
    value = 0.0;
    double stabw = 0.0;
    int i = 0;

    // positions
    positions R1_new =positions(parameters);
    positions * Rnew = &R1_new;

    positions R1_old =positions(parameters);
    positions * Rold = &R1_old;

    double Phi = fct->getValue(Rnew);
    double P_old = Phi*Phi;
    double P_new;


    while (i<nSamples)
    {


        // perform step
        mat newPositions = Rold->get_pos()+ (+2*randu<mat>(ndim,nParticles)-ones(ndim,nParticles))*stepSize;
        Rnew->set_pos(newPositions);

        // calculate P_new/P_old
         Phi = fct->getValue(Rnew);
         P_new =Phi*Phi;



        // compare to random variable
         vec zufallszahl(randu(1));
         if( zufallszahl(0) <= P_new/P_old )
        {
            P_old = P_new;
            swap(Rnew,Rold);

            // increase accepted steps
            acceptedSteps++;
        }
        else
        {
            // keep old position
        }

        // add energy value
         double ede = H->localEnergy(fct,Rold);
         value += ede;
         stabw += ede*ede;

        // increase i
        i++;

    }

    value /= nSamples*nParticles;
    variance = sqrt(stabw-value*value) / (nParticles*nSamples);
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
    return (double) acceptedSteps/nSamples;
}

double mcInt::get_value()
{
    return value;
}


double mcInt::get_variance()
{
    return variance;
}
