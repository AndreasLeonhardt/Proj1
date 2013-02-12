#include <iostream>
#include <armadillo>
#include <libconfig.h++>
#include <string>

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

    ofstream results;
    results.open("Results_Proj1.txt");
    results << "alpha\tbeta\tE\tdE\tE_analytical\tdE_analytical\tacceptance_ratio" << endl ;

    //loop over different parameters alpha, beta
    double a_min = 3.0;
    double a_sup = 3.7;
    double a_step= 0.2;
    int a_n = floor((a_sup-a_min)/a_step)+1;

    double b_min = .5;
    double b_sup = 1.6;
    double b_step= .5;
    int b_n = floor((b_sup-b_min)/b_step)+1;
    int c_n=b_n*a_n;
    int c = 0;

    for (double a = a_min; a<a_sup ; a+=a_step)
    {
        for(double b = b_min; b<b_sup; b+=b_step)

        {
            fun->set_alpha(a);
            fun->set_beta(b);

            MC.integrate(fun,H,parameters);

            results << a << "\t"
                    << b << "\t"
                    << MC.get_value()  << "\t"
                    << MC.get_variance() << "\t"
                    << MC.get_value_analytical() << "\t"
                    << MC.get_variance_analytical() << "\t"
                    << MC.get_acceptanceRatio()*100 << "%"
                    <<endl;


             c++;


             int counter = 100*c/c_n;

             for (int i=0;i<45;i++)
                 cout<<"\n\n\n\n"<<endl;

             string bar(counter/1.5151,'X');
             bar.append((100-counter)/1.5151,'_');

             cout << "progress = " << counter << "%" << endl;
             cout << bar << endl;

        }

    }

    results.close();

    return 0;
}




// EOF
