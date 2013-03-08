#include <iostream>
#include <armadillo>
#include <libconfig.h++>
#include <string>

#include "positions.h"
#include "function.h"
#include "trialfct.h"
#include "trialfct_analytical.h"
#include "hamilton.h"
#include "hamilton_analytical.h"
#include "hamilton_numerical.h"
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


    // seed for rand0
    long int idum = -1;


    // create instance of TrialFct
    function * fun = new function( parameters);


    // create instance of hamilton
    hamilton * H = new hamilton(parameters);

    int schalter = parameters->lookup("analytical_energy_density");
    if (!schalter)
    {   
        H = new hamilton_numerical(parameters);
        fun = new TrialFct(parameters);
    }
    else
    {
        H = new hamilton_analytical(parameters);
        //fun = new TrialFct_analytical(parameters);
        fun = new TrialFct(parameters);

    }


    // create instance of mcInt
    mcInt MC = mcInt(parameters);

    ofstream results;
    results.open(parameters->lookup("outputfile"));
    results << "Integration points: " << (int) parameters->lookup("nSamples") << "  Z=" << (int) parameters->lookup("Z")
            << "  analytical: " << (int) parameters->lookup("analytical_energy_density") <<endl;
    results << "alpha\tbeta\tE\tdE\tacceptance_ratio" << endl ;

    //loop over different parameters alpha, beta
    double a_min = parameters->lookup("a_min");
    double a_step = parameters->lookup("a_step");
    double a_max = parameters->lookup("a_max");
    a_max += 0.000001*a_step; // rather dirty way to get rid of rounding errors.
    double a_n = floor((a_max-a_min)/a_step) +1;

    double b_min = parameters->lookup("b_min");
    double b_step = parameters->lookup("b_step");
    double b_max = parameters->lookup("b_max");
    b_max += 0.000001*b_step; // rather dirty way to get rid of rounding errors.
    int b_n = floor((b_max-b_min)/b_step) +1;

    double c_n=b_n*a_n;
    int c = 0;
    string bar;

    for (double a = a_min; a<a_max ; a+=a_step)
    {
        for(double b = b_min; b<b_max; b+=b_step)

        {
            fun->setParameter(a,0);
            fun->setParameter(b,1);
            positions * Rinitial = MC.thermalise(fun,idum, parameters);
            MC.integrate(fun,H, Rinitial,idum, parameters);

            results << a << "\t"
                    << b << "\t"
                    << MC.get_value()  << "\t"
                    << MC.get_variance() << "\t"
                    << MC.get_acceptanceRatio()*100 << "%"
                    <<endl;

            // status of calculation
            c++;
            int counter = 100*c/c_n;
            // clear screen (kind of)
            for (int i=0;i<10;i++)
                 cout<<"\n\n\n\n"<<endl;
            // create string for status bar
            bar.assign(counter/1.5151,'X');
            bar.append((100-counter)/1.5151,'_');
            // write status, percentage and status bar
            cout << "progress = " << counter << "%" << endl;
            cout << bar << endl;

        }

    }

    delete fun;
    delete H;
    results.close();

    return 0;
}




// EOF
