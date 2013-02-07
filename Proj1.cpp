#include <iostream>
#include <armadillo>
#include <libconfig.h++>

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

    // the position of particle i is stored in the i-th column.
    // get matrix size
    int ndim = parameters->lookup("ndim");
    int nParticles = parameters->lookup("nParticles");
    // create position matrix with random numbers
    mat positions =randu(ndim,nParticles);

    // parameters of trial wave function
    double alpha = parameters->lookup("alpha");
    double beta = parameters->lookup("beta");


    // create instance of TrialFct
    TrialFct functio = TrialFct( parameters, alpha, beta );
    TrialFct * fun = &functio;
    fun->set_position(positions);

    //
    hamilton  Ham = hamilton(parameters);
    hamilton * H = &Ham;
    // create instance of mcInt
    mcInt MC = mcInt(parameters);
    MC.integrate(fun,H);



    cout << "E=" << MC.get_value() << endl <<"  dE=" << MC.get_variance() << endl << "acceptance ratio=" << MC.get_acceptanceRatio()*100<< "%" <<endl;

    return 0;
}




// EOF
