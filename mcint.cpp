#include "mcint.h"

mcInt::mcInt()
{
    nSamples = 1000;
    timeStep = 1;
    ndim = 3;
    nParticles = 2;
    sqrtTimeStep = 1;
}


mcInt::mcInt(Config * parameters)
{
    nSamples = parameters->lookup("nSamples");
    timeStep = parameters->lookup("timeStep");
    ndim = parameters->lookup("ndim");
    nParticles = parameters->lookup("nParticles");
    sqrtTimeStep = sqrt(timeStep);
}




void mcInt::integrate(function * fct, hamilton * H, Config * parameters, long int idum)
{
    acceptedSteps = 0;
    value = 0.0;
    double stabw = 0.0;
    int n = 0;

    // positions
    positions R1_new =positions(parameters);
    positions * Rnew = &R1_new;

    positions R1_old =positions(parameters);
    positions * Rold = &R1_old;

    double Phi = fct->getValue(Rnew);
    double P_old = Phi*Phi;
    double P_new;




    while (n<nSamples)
    {


        // perform step one particles at the time
        for (int i =0; i<nParticles; i++)
        {
            vec newPosition(ndim);
            newPosition = Rold->get_r(i)+randn(ndim)*sqrtTimeStep +0.5*fct->quantumForce(i,Rold)*timeStep;
            Rnew->set_singlePos(newPosition,i);

            // calculate Greensfct CHECK FOR SIGN ERROR
            vec Fold = fct->quantumForce(i,Rold);
            vec Fnew = fct->quantumForce(i,Rnew);

            double ratioGreensfunction = 0.5*dot(
                                            (Fold+Fnew),
                                            ( Rold->get_r(i) - Rnew->get_r(i) + 0.5*timeStep*(Fold-Fnew)) // there might be a sign errror
                                            );
            ratioGreensfunction = exp(-ratioGreensfunction);


            // calculate ( G(old,new)*P_new ) / ( G(new,old)*P_old )
            Phi = fct->getValue(Rnew);
            P_new =Phi*Phi;


            double Zufallszahl = ran0(&idum);

            if( Zufallszahl <= ratioGreensfunction*P_new/P_old )
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
        } // end loop over particles.


        // add energy value
         double ede = H->localEnergy(fct,Rold);
         value += ede;
         stabw += ede*ede;



        // increase counter
        n++;

    }

    value /= nSamples;
    variance = sqrt(stabw-value*value) / nSamples;
    acceptedSteps/=nParticles;
}




void mcInt::set_nSamples(int NewnSamples)
{
    nSamples = NewnSamples;
}

int  mcInt::get_nSamples()
{
    return nSamples;
}

void mcInt::set_timeStep(double NewTimeStep)
{
    timeStep=NewTimeStep;
}

int  mcInt::get_timeStep()
{
    return timeStep;
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

