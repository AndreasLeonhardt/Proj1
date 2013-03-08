#include "mcint.h"

mcInt::mcInt()
{
    nSamples = 1000;
    thermalisationSteps = 100;
    timeStep = 1;
    ndim = 3;
    nParticles = 2;
    sqrtTimeStep = 1;
}


mcInt::mcInt(Config * parameters)
{
    nSamples = parameters->lookup("nSamples");
    thermalisationSteps = parameters->lookup("thermalisationSteps");
    timeStep = parameters->lookup("timeStep");
    ndim = parameters->lookup("ndim");
    nParticles = parameters->lookup("nParticles");
    sqrtTimeStep = sqrt(timeStep);
}





positions * mcInt::Step(function * fct,  positions * Rold, long int idum, Config * parameters)
{
    double P_new;
    positions * Rnew = new positions(parameters);
    vec newPosition(ndim);

    Rnew->set_pos(Rold->get_pos());

    // perform step one particles at the time
    for (int i =0; i<nParticles; i++)
    {
        newPosition = Rold->get_singlePos(i)+randu(ndim)*sqrtTimeStep +0*0.5*fct->quantumForce(i,Rold)*timeStep;

        Rnew->set_singlePos(newPosition,i);

        // calculate Greensfct CHECK FOR SIGN ERROR
        vec Fold(ndim);
        Fold = fct->quantumForce(i,Rold);
        vec Fnew(ndim);
        Fnew = fct->quantumForce(i,Rnew);

        // there might be a sign errror in calculating the exponent
        double ratioGreensfunction = 0.5*dot(
                                             (Fold+Fnew),
                                             ( Rold->get_singlePos(i) - Rnew->get_singlePos(i) + 0.5*timeStep*(Fold-Fnew))
                                            );

        ratioGreensfunction = exp(ratioGreensfunction);
        ratioGreensfunction =1;

        // calculate ( G(old,new)*P_new ) / ( G(new,old)*P_old )
        double Phi = fct->getValue(Rnew);
        P_new =Phi*Phi;
// cout << P_new << endl;

        double Zufallszahl = ran0(&idum);

        if( Zufallszahl <= ratioGreensfunction*P_new/P_old )
        {
            P_old = P_new;
           Rold->set_singlePos(newPosition,i);
            //Rold->set_pos(Rnew->get_pos())

            // increase accepted steps
            acceptedSteps++;
        }
        else
        {
            Rnew->set_singlePos(Rold->get_singlePos(i),i);
           // Rnew->set_pos(Rold->get_pos());
        }


    } // end loop over particles.

    return Rold;

}


positions * mcInt::thermalise(function * fct, long int idum, Config * parameters)
{
    acceptedSteps=0;

    positions * Rold = new positions(parameters);

    P_old = fct->getValue(Rold)*fct->getValue(Rold);

    for (int i=0;i<thermalisationSteps;i++)
    {
        Rold=Step(fct, Rold,idum, parameters);
    }

    return Rold;
}


void mcInt::integrate(function * fct, hamilton * H, positions *Rold, long int idum, Config * parameters)
{

    value = 0.0;
    double stabw = 0.0;

    double Phi = fct->getValue(Rold);
     P_old = Phi*Phi;
    for(int n = 0; n<nSamples; n++)
    {
        // perform step
        Rold = Step(fct, Rold, idum, parameters);


        // add energy value
         double ede = H->localEnergy(fct,Rold);

         // test for hydrogen case
         double rrr = Rold->get_r(0);
         ede = -1/rrr - fct->getParameter(0)*(0.5*fct->getParameter(0)-1/rrr);

         value += ede;
         stabw += ede*ede;

        // cout << ede  << " r = " << rrr << "   parameter: " << fct->getParameter(0) << endl;
    }

    value /= nSamples;
    stabw /= nSamples;

    variance = sqrt ((stabw-value*value) / nSamples);


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
    return (double) acceptedSteps/(nSamples+thermalisationSteps);
}

double mcInt::get_value()
{
    return value;
}


double mcInt::get_variance()
{
    return variance;
}

