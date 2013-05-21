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





// performs the Monte Carlo integration
// calculating the mean and writing the samples to a file for later analysis.
// direct calculation might be removed due to speed issue and inaccuracy concerning the error.
void mcInt::integrate(function * fct, hamilton * H, long * idumadress,Config * parameters,char * samplefile)
{
    // opening file for sample storage
    ofstream outfile;


    outfile.open(samplefile,ios::out | ios::binary);

    acceptedSteps=0;

    positions * R = new positions(parameters);

    fct->setSlaterinv(R);


    for (int i=0;i<thermalisationSteps;i++)
    {
        R=Step(fct, R,idumadress,parameters);
    }

    // store samples to write blockwise
    double samples[10000];
    // loop over steps in blocks of 10000
    for(int n=0;n<nSamples/10000;n++)
    {
        for(int m = 0; m<10000; m++)
        {
        // perform step
        R = Step(fct, R, idumadress,parameters);

        // add energy value
         double ede = H->localEnergy(fct,R);
         samples[m]=ede;
        }

        // write sample to file
        outfile.write((char*) &samples, sizeof(double)*10000);
    }
    // do the rest, that didn't fit in the last block
    for (int n=0;n<nSamples%10000;n++)
    {
        // perform step
        R = Step(fct, R, idumadress,parameters);

        // add energy value
         double ede = H->localEnergy(fct,R);
         samples[n]=ede;

    }


    // write sample to file but only the new values
    outfile.write((char*) &samples, sizeof(double)*nSamples%10000);

    acceptedSteps/=nParticles;

    // close sample close
    outfile.close();
}




// function StatGrad (statistical gradient)
// returns the gradient of the trial wave function according to the parameters
// R0 is not treated as a parameter any longer, because it is a part of the system
// and not only the trial wave funciton.
// the 'statistical' is due to the big noise in the Monte Carlo calculation
// because only a few cycles are used. This is actually used in the statistical gradient approach.
vec mcInt::StatGrad(function * fct, hamilton *H,long * idumadress,int nParams, Config *parameters)
{


    vec store1 = zeros(nParams);
    vec store2 = zeros(nParams);
    double storeEnergy=0.0;
    double E,derfct;


    // thermalization  (could avoided if old position is used again for the next optimization step.)
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
    vec result =2.0/SCnSamples*(store2-store1*storeEnergy);

    return result;

}




mat mcInt::blocking(Config *parameters,char * samplefilebody, int rmax)
{

    int bmin = parameters->lookup("minBlockSize");
    int bmax = parameters->lookup("maxBlockSize");
    int bsteps= parameters->lookup("BlockSteps");



    // opening file for sample storage
    char * file = new char[50];
    ifstream infile;
    int NumberOfSamples=0;

    // when running program in paralell, we need different names here.
    // (adding numbers and so on)
            struct stat fileproperties;
    for (int r=0;r<rmax;r++)
    {
        sprintf(file,"%s_%u.dat",samplefilebody,r);
        // get size of data block


        if (stat(file,&fileproperties)==0)
        {
            NumberOfSamples += fileproperties.st_size/sizeof(double);
        }

    }

    // allocate data block
    double data[NumberOfSamples];

    for (int r=0;r<rmax;r++)
    {
        sprintf(file,"%s_%u.dat",samplefilebody,r);
        stat(file,&fileproperties);
        infile.open(file,ios::in | ios::binary);
        // write data into array
        infile.read((char*)&data,fileproperties.st_size);
        infile.close();
    }

    value = 0.0;
    for(int i=0;i<NumberOfSamples;i++)
    {
        value+=data[i];
    }
    value /=NumberOfSamples;


    // intialize useful stuff for the loop
    mat std=zeros(2,bsteps);
    double newSample;
    double av,sqrav;
    int b=bmin;
    double bb=bmin;
    int n;

    // loop over different block sizes
    double blockstepsize = (double) (bmax-bmin)/bsteps;

    for (int i=0;i<bsteps;i++)
    {
        // get number of blocks
        n= NumberOfSamples/b;
        // reset storage for average and squareaverage
        av=0.0;
        sqrav=0.0;

        // loop over all blocks
        for (int j=0;j<n;j++)
        {
            newSample=0.0;
            // calculate averages over blocks
            for(int i = 0;i<b;i++)
            {
                newSample+=data[j*b+i];
            }
            newSample/=b;

            // add value of each block
            av+=newSample;
            sqrav+=newSample*newSample;
        }

        av/=n;
        sqrav/=n;
        std(0,i)=b;
        std(1,i)=sqrt( (sqrav-av*av)/n );
        bb+=blockstepsize;
        b=round(bb);
    }

    return std;

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

