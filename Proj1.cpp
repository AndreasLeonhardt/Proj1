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
    int parameterIterations = parameters->lookup("parameterIterations");

        vec a=zeros(nParams);
        for (int i=0;i<nParams;i++)
        {
            a(i)=parameters->lookup("Parameters.[i]");
            fun->setParameter(a(i),i);
        }
        // set alpha to Z for a start
        int Z = parameters->lookup("Z");
        a(0)=Z;
        fun->setParameter(a(0),0);


        // without adaptive stepsize, using 1/i
        for(int i=1;i<parameterIterations+1;i++)
        {
            a -= MC.StatGrad(fun,H,idumadress,nParams,parameters)/i;
            fun->setParameter(a);
            cout<< "Parameters: "<<a<<endl;
            MC.integrate(fun,H,idumadress,parameters);
            results << a(0)<<"\t"<<a(1)<<"\t"<<MC.get_value()<<"\t"<<MC.get_variance()<<"\t"<<MC.get_acceptanceRatio()<<endl;
            cout <<"energy: "<<MC.get_value()<<endl;

        }





    delete fun;
    delete H;
    results.close();

    return 0;
}




// EOF
