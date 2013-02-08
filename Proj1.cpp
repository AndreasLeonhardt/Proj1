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

    //
    hamilton  Ham = hamilton(parameters);
    hamilton * H = &Ham;
    // create instance of mcInt
    mcInt MC = mcInt(parameters);

    MC.integrate(fun,H,parameters);


    cout << "hallo"<< endl;
    cout << "E=" << MC.get_value() << endl <<"  dE=" << MC.get_variance() << endl << "acceptance ratio=" << MC.get_acceptanceRatio()*100<< "%" <<endl;

    return 0;
}




// EOF
