#include "mcint.h"

mcInt::mcInt()
{
    nSamples = 1000;
    thermalisationSteps = 100;
    SCnSamples=1000;
    SCthermalisationSteps =100;
    timeStep = 1;
    ndim = 3;
    nParticles = 2;
    sqrtTimeStep = 1;
}


mcInt::mcInt(Config * parameters)
{
    nSamples = parameters->lookup("nSamples");
    thermalisationSteps = parameters->lookup("thermalisationSteps");
    SCnSamples=parameters->lookup("SCnSamples");
    SCthermalisationSteps = parameters->lookup("SCthermalisationSteps");
    timeStep = parameters->lookup("timeStep");
    ndim = parameters->lookup("ndim");
    nParticles = parameters->lookup("nParticles");
    sqrtTimeStep = sqrt(timeStep);
}





positions * mcInt::Step(function * fct,  positions * Rold, long int * idumadress, Config * parameters)
{
    positions Rnewposition = positions(*Rold);
    positions * Rnew = &Rnewposition;
    vec newPosition(ndim);
    vec Fold(ndim);
    vec Fnew(ndim);

    // perform step one particles at the time
    for (int i =0; i<nParticles; i++)
    {
        // calculate quantum Force at old position
        Fold = fct->quantumForce(i,Rold);


        // calculate proposal for new position
        newPosition = Rold->get_singlePos(i) + randn(ndim)*sqrtTimeStep +0.5*Fold*timeStep;
        Rnew->set_singlePos(newPosition,i);

        // calculate the shortcut for the ratio
        double ratioSlater = fct->SlaterRatio(i,Rold,Rnew);
        double ratio = ratioSlater*fct->JastrowRatio(i,Rold,Rnew);

        // update inverse Slater matrix.
        // save old Slaterinv first, in case the step is not accepted.
        mat oldSlaterinv = fct->getinvslatermatrix(i);
        fct->updateSlaterinv(i,Rnew,ratioSlater);
        // calculate quantum force at new position,
        Fnew = fct->quantumForce(i,Rnew);


        // Greens function CHECK FOR SIGN ERROR
        // there might be a sign errror in calculating the exponent
        double ratioGreensfunction = 0.5*dot(
                                             (Fold+Fnew),
                                             ( Rold->get_singlePos(i) - Rnew->get_singlePos(i) + 0.5*timeStep*(Fold-Fnew))
                                            );

        ratioGreensfunction = exp(ratioGreensfunction);




        // test acceptance
        double Zufallszahl = ran0(idumadress);
        if( Zufallszahl <= ratioGreensfunction*ratio*ratio )
        {

           Rold->set_singlePos(newPosition,i);
            // increase accepted steps
            acceptedSteps++;
        }
        else
        {
            Rnew->set_singlePos(Rold->get_singlePos(i),i);
            fct->setSlaterinv(i,oldSlaterinv);
        }

    } // end loop over particles.

    return Rold;
}






void mcInt::integrate(function * fct, hamilton * H, long * idumadress,Config * parameters)
{


    acceptedSteps=0;

    positions * R = new positions(parameters);

    fct->setSlaterinv(R);


    for (int i=0;i<thermalisationSteps;i++)
    {
        R=Step(fct, R,idumadress,parameters);
    }

    value = 0.0;
    double stabw = 0.0;

    for(int n = 0; n<nSamples; n++)
    {
        // perform step
        R = Step(fct, R, idumadress,parameters);

        // add energy value
         double ede = H->localEnergy(fct,R);




         value += ede;
         stabw += ede*ede;
    }

    value /= nSamples;
    stabw /= nSamples;

    variance = sqrt ((stabw-value*value) / nSamples);


    acceptedSteps/=nParticles;
}





vec mcInt::StatGrad(function * fct, hamilton *H,long * idumadress,int nParams, Config *parameters)
{


    vec store1 = zeros(nParams);
    vec store2 = zeros(nParams);
    double storeEnergy=0.0;
    double E,derfct;


    // thermalization
    positions * R = new positions(parameters);
    fct->setSlaterinv(R);
    for (int n=0;n<SCthermalisationSteps;n++)
    {
        R=Step(fct, R,idumadress,parameters);
    }


    // mini MC integration
    for(int n = 0; n<SCnSamples; n++)
    {
        // perform step
        R = Step(fct, R, idumadress,parameters);

        // add energy value
         E = H->localEnergy(fct,R);
         for (int d=0;d<nParams;d++)
         {
             derfct=fct->ParamDerivativeOverFct(R,d);
             store1[d]+=derfct;
             store2[d]+=derfct*E;
         }

         storeEnergy +=E;
    }

    storeEnergy/=SCnSamples;
    return 2.0/SCnSamples*(store2-store1*storeEnergy);
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

