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


    // create instance of TrialFct
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
    results << "Integration points: " << (int) parameters->lookup("nSamples") << "  Z=" << (int) parameters->lookup("Z")
            << "  analytical: " << (int) parameters->lookup("analytical_energy_density") <<endl;
    results << "alpha\tbeta\tE\tdE\tacceptance_ratio" << endl ;

    //loop over different parameters alpha, bet
    int nParams = parameters->lookup("nParameters");
//    double resLimit = parameters->lookup("ResidualLimit");
//    double NewtonLimit = parameters->lookup("NewtonLimit");
//    int Newtoncounterlimit =parameters->lookup("NewtonCounterLimit");

    vec a=zeros(nParams);
        vec gradient(nParams);
        for (int i=0;i<nParams;i++)
        {
            a(i)=parameters->lookup("Parameters.[i]");
        }


        // without adaptive stepsize
        for(int i=1;i<parameterIterations;i++)
        {

            a -= MC.StatGrad(fun,H,idumadress,parameters)/i;
        }



    delete fun;
    delete H;
    results.close();

    return 0;
}




// EOF
