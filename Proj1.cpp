#include <iostream>
#include <armadillo>
#include <libconfig.h++>
#include <string>

#include "positions.h"
#include "function.h"
#include "trialfct.h"
#include "trialfct_analytical.h"
#include "hamilton.h"
#include "mcint.h"

using namespace std;
using namespace arma;
using namespace libconfig;





int main()
{
    // INITIALIZATION ----------------------------------------------------------------------
    cout << "initialization "<<flush;

    // read parameters from file
    Config conf_parameters;
    Config * parameters = &conf_parameters;
    parameters->readFile("../Proj1/parameters.cfg");

    int checkforevenparticlenumber = parameters->lookup("nParticles");
    checkforevenparticlenumber = checkforevenparticlenumber%2;
    if (checkforevenparticlenumber)
    {
        cout << "Programm might not run properly for an odd number of Particles." << endl
                << "Change  \"nParticles\" in \"parameters.config\" to an even number" << endl;
        cout << checkforevenparticlenumber <<endl;
    }

    // seed for rand0
    long  idum;
    idum =-92;
    long int * idumadress = &idum;


    // create pointer to function
     function * fun;
    // create instance of hamilton
     hamilton * H =new hamilton(parameters);

    int schalter = parameters->lookup("analytical_energy_density");
    if (!schalter)
    {   
        fun = new TrialFct(parameters);

    }
    else
    {
        fun = new TrialFct_analytical(parameters);
    }

    // create instance of mcInt
    mcInt MC = mcInt(parameters);


    ofstream results;
    results.open(parameters->lookup("outputfile"));
    results << "Integration points: " << (int) parameters->lookup("nSamples")
            << "  analytical: " << (int) parameters->lookup("analytical_energy_density") <<endl;
    results << "alpha\tbeta\tR_0\tE\tacceptance_ratio" << endl ;

    //loop over different parameters alpha, beta and R_0
    int nParams = parameters->lookup("nParameters");
    int parameterIterations = parameters->lookup("parameterIterations");

        vec a=zeros(nParams);
        a(0)=parameters->lookup("Parameters.[0]");
        a(1)=parameters->lookup("Parameters.[1]");

        cout << "..... done."<<endl;


        // PARAMETER OPTIMIZATION ---------------------------------------------------------------------------
        cout<<"parameter optimization "<<flush;

        // set alpha to Z for a start
        int Z = parameters->lookup("Z.[0]");
        a(0)=Z;
        fun->setParameter(a(0),0);



        // without adaptive stepsize, using 1/i as factor
        for(int i=1;i<parameterIterations+1;i++)
        {
            a -= MC.StatGrad(fun,H,idumadress,nParams,parameters)/i;          
                fun->setParameter(a);
        }

        cout << "..... done."<<endl;


        // INTEGRATION ---------------------------------------------------------------------------------------
        cout <<"integration "<<endl;

        // perform Monte Carlo integration
        MC.integrate(fun,H,idumadress,parameters);


        // BLOCKING AND WRITING OF RESULTS --------------------------------------------------------------------
        // analyse the result via blocking
        cout<<"blocking "<<endl;
        mat blockingResult = MC.blocking(parameters).t();



        // write results
        cout<<"write results "<<endl;
        results << a(0)<<"\t"<<a(1)<<"\t"<<fun->get_R0()<<"\t"<<MC.get_value()<<"\t"<<"\t"<<MC.get_acceptanceRatio()<<endl<<endl;

        results<<"blocking result:"<<endl
               <<"blocksize\tvariance"<<endl
               <<blockingResult<<endl;


    delete fun;
    delete H;
    results.close();

    // might call a (python)script from here that plots the data in the result file
    // using matplotlib for python for example
    // call the script using the system() function

    char* command=new char[50];
    const char* outputfile = parameters->lookup("outputfile");
    sprintf(command,"python ../Proj1/plot.py %s &", outputfile);

    system(command);

    return 0;
}




// EOF
