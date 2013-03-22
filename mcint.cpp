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





positions * mcInt::Step(function * fct,  positions * Rold, long int * idumadress)
{
    positions Rnewposition = positions(*Rold);
    positions * Rnew = &Rnewposition;
    vec newPosition(ndim);
    vec Fold(ndim);
    vec Fnew(ndim);


// compare two functions
    TrialFct* fctII = new TrialFct();
    fctII->setParameter(fct->getParameter(0),0);
    fctII->setParameter(fct->getParameter(1),1);
    fctII->setnParticles(fct->getnParticles());



    // perform step one particles at the time
    for (int i =0; i<nParticles; i++)
    {
        // calculate quantum Force at old position
        Fold = fct->quantumForce(i,Rold);
        // test difference
        //cout << "quantum force analytical: "<<endl<<Fold<<"qf  num: "<<endl<<fctII->quantumForce(i,Rold)<<endl;
        //cout << i<<"================" <<endl;
        // calculate proposal for new position
        newPosition = Rold->get_singlePos(i) + randn(ndim)*sqrtTimeStep +0.5*Fold*timeStep;
        Rnew->set_singlePos(newPosition,i);

        // calculate quantum force at new position
        Fnew = fct->quantumForce(i,Rnew);
        //cout << "quantum force: "<<endl<< Fold <<endl <<Fnew<<endl;

        // Greens function CHECK FOR SIGN ERROR
        // there might be a sign errror in calculating the exponent
        double ratioGreensfunction = -0.5*dot(
                                             (Fold+Fnew),
                                             ( Rold->get_singlePos(i) - Rnew->get_singlePos(i) + 0.5*timeStep*(Fold-Fnew))
                                            );

        ratioGreensfunction = exp(ratioGreensfunction);



        // NEW: calculate the shortcut for the ratio
        double ratioSlater = fct->SlaterRatio(i,Rold,Rnew);
        //double ratioJastrow = fct->JastrowRatio(i,Rold,Rnew); without Jastrow factor for a start

        // test acceptance
        double Zufallszahl = ran0(idumadress);
        if( Zufallszahl <= ratioGreensfunction*ratioSlater*ratioSlater )
        {

           Rold->set_singlePos(newPosition,i);


           // update inverse slater matrix
           //fct->updateSlaterinv(i,Rold,ratioSlater);
           //mat test = fct->getinvslatermatrix(1);
           fct->setSlaterinv(Rold);
           //mat testII = fctII->getinvslatermatrix(1);
           //cout << test/testII<<endl;

            // increase accepted steps
            acceptedSteps++;
        }
        else
        {
            Rnew->set_singlePos(Rold->get_singlePos(i),i);
           // Rnew->set_pos(Rold->get_pos());
        }


    } // end loop over particles.
    delete fctII;
    return Rold;


}


positions * mcInt::thermalise(function * fct, long int * idumadress, Config * parameters)
{
    acceptedSteps=0;

    positions * Rold = new positions(parameters);

    fct->setSlaterinv(Rold);


    for (int i=0;i<thermalisationSteps;i++)
    {
        Rold=Step(fct, Rold,idumadress);
    }

    return Rold;
}


void mcInt::integrate(function * fct, hamilton * H, positions *Rold, long * idumadress)
{

    value = 0.0;
    double stabw = 0.0;

    for(int n = 0; n<nSamples; n++)
    {
        // perform step
        Rold = Step(fct, Rold, idumadress);


        // add energy value
         double ede = H->localEnergy(fct,Rold);

         value += ede;
         stabw += ede*ede;
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

