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





positions * mcInt::Step(function * fct,  positions * Rold, long int * idumadress, Config * parameters)
{
    double P_new;
    positions * Rnew = new positions(parameters, idumadress);
    vec newPosition(ndim);


    Rnew->set_pos(Rold->get_pos());

    // perform step one particles at the time
    for (int i =0; i<nParticles; i++)
    {
        // calculate quantum Force at old position
        vec Fold(ndim);
        Fold = fct->quantumForce(i,Rold);

        // calculate proposal for new position
        newPosition = Rold->get_singlePos(i) + randn(ndim)*sqrtTimeStep +0.5*Fold*timeStep;
        Rnew->set_singlePos(newPosition,i);

        // calculate quantum force at new position
        vec Fnew(ndim);
        Fnew = fct->quantumForce(i,Rnew);

        // Greens function CHECK FOR SIGN ERROR
        // there might be a sign errror in calculating the exponent
        double ratioGreensfunction = 0.5*dot(
                                             (Fold+Fnew),
                                             ( Rold->get_singlePos(i) - Rnew->get_singlePos(i) + 0.5*timeStep*(Fold-Fnew))
                                            );

        ratioGreensfunction = exp(ratioGreensfunction);

        // calculate probability density at the new position
        double Phi = fct->getValue(Rnew);
        P_new =Phi*Phi;

        // test acceptance
        double Zufallszahl = ran0(idumadress);
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

/*
       double P_new;
        positions * Rnew = new positions(parameters, idumadress);
        mat randompart=zeros(ndim,nParticles);

        for (int j=0;j<nParticles;j++)
        {
            for (int i = 0; i<ndim; i++)
            {   double number = ran0(idumadress);
                randompart(i,j)=number-0.5;
                //cout << number << endl;
            }
        }


        Rnew->set_pos(Rold->get_pos()+randompart);


            double a= fct->getParameter(0);
            double rnew= Rnew->get_r(0);
            double rold = Rold->get_r(0);
            // calculate ( G(old,new)*P_new ) / ( G(new,old)*P_old )
            double Phi = fct->getValue(Rnew);
            // test with closed form expression for hydrogen case:
            Phi = exp(-a*rnew);
            P_new =Phi*Phi;
            double Phiold = exp(-a*rold);
            P_old = Phiold*Phiold;




            double Zufallszahl = ran0(idumadress);

            if( Zufallszahl <= P_new/P_old )
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



        return Rold;


*/
}


positions * mcInt::thermalise(function * fct, long int * idumadress, Config * parameters)
{
    acceptedSteps=0;

    positions * Rold = new positions(parameters, idumadress);

    P_old = fct->getValue(Rold)*fct->getValue(Rold);

    for (int i=0;i<thermalisationSteps;i++)
    {
        Rold=Step(fct, Rold,idumadress, parameters);
    }

    return Rold;
}


void mcInt::integrate(function * fct, hamilton * H, positions *Rold, long * idumadress, Config * parameters)
{

    value = 0.0;
    double stabw = 0.0;

    double Phi = fct->getValue(Rold);
     P_old = Phi*Phi;
    for(int n = 0; n<nSamples; n++)
    {
        // perform step
        Rold = Step(fct, Rold, idumadress, parameters);


        // add energy value

        // THE VALUE FROM THE LOCAL ENERGY SEEMS TO BE WRONG.
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

