#include <iostream>
#include <armadillo>
#include <libconfig.h++>

#include "positions.h"
#include "trialfct.h"
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
    parameters->readFile("parameters.cfg");

    // seed for randu
    srand(41);

    // create instance of TrialFct
    TrialFct functio = TrialFct( parameters);
    TrialFct * fun = &functio;

    // create instance of hamilton
    hamilton  Ham = hamilton(parameters);
    hamilton * H = &Ham;

    // create instance of mcInt
    mcInt MC = mcInt(parameters);

    //loop over different parameters alpha, beta
    for (double a = 2.5; a<4.6 ; a+=0.5)
    {
        for(double b = 1.0; b<4.1; b+=1.0)
        {
            fun->set_alpha(a);
            fun->set_beta(b);

            MC.integrate(fun,H,parameters);

            cout << "alpha = " << a << "  beta = " << b << endl;
            cout << "E=" << MC.get_value() << endl <<"dE=" << MC.get_variance() << endl;
            cout << "acceptance ratio=" << MC.get_acceptanceRatio()*100<< "%" <<endl << endl;

        }
    }



    return 0;
}




// EOF
