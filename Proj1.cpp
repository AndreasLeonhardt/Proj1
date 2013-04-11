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
    double resLimit = parameters->lookup("ResidualLimit");
    double NewtonLimit = parameters->lookup("NewtonLimit");
    int Newtoncounterlimit =parameters->lookup("NewtonCounterLimit");

        double a[nParams];
        double here;
        vec forward(nParams);
        double E,Em,Ep;
        vec d(nParams);
        vec r(nParams);
        vec gradient(nParams);
        double s,sold, t;
        double resLength =1000;
        double Newtondiff = 1000;
        for (int i=0;i<nParams;i++)
        {
            a[i]=parameters->lookup("Parameters.[i]");
        }
    int maxcounter = parameters->lookup("maxIterations");
    int counter = 0;
    int Newtoncounter;
    double h = parameters->lookup("paramstepwidth");


    // set initial conditions
    // calculate -gradient
    MC.integrate(fun,H,MC.thermalise(fun,idumadress,parameters),idumadress,parameters);
    here = MC.get_value();
    for(int i=0;i<nParams;i++)
    {
        // move on step forward
        fun->setParameter(a[i]+h,i);
        MC.integrate(fun,H,MC.thermalise(fun,idumadress,parameters),idumadress,parameters);
        forward[i]=MC.get_value();

        // set back to initial value
        fun->setParameter(a[i],i);
    }



    r = ( -forward + ones(nParams)*here )/h;
    d=r;

    while (counter<maxcounter)
    {

        cout<<"outer loop, cycle number: "<<counter<<endl;

        // set parameters
        for(int i=0;i<nParams;i++)
        {
            fun->setParameter(a[i],i);
        }

        // find initial position thorugh thermalisation
        positions * Rinitial = MC.thermalise(fun, idumadress, parameters);
        // actual calculation
        MC.integrate(fun,H, Rinitial,idumadress,parameters);

        // write results
        results << a[1] << "\t"
                << a[2] << "\t"
                << MC.get_value()  << "\t"
                << MC.get_variance() << "\t"
                << MC.get_acceptanceRatio()*100 << "%"
                <<endl;


        // minimize E(x+s*r)
        s=0.0;
        sold = 1000.0;
        Newtoncounter=0;

        while (Newtondiff>NewtonLimit && Newtoncounter<Newtoncounterlimit)
        {
            cout<<"Inner loop, cycle number: "<<Newtoncounter<<endl;

            for(int i=0;i<nParams;i++)
            {

                MC.integrate(fun,H,MC.thermalise(fun,idumadress,parameters),idumadress,parameters);
                E=MC.get_value();

                fun->setParameter(a[i]+h*d[i],i);
                MC.integrate(fun,H,MC.thermalise(fun,idumadress,parameters),idumadress,parameters);
                Ep = MC.get_value();

                fun->setParameter(a[i]-h*d[i],i);
                MC.integrate(fun,H,MC.thermalise(fun,idumadress,parameters),idumadress,parameters);
                Em = MC.get_value();

                fun->setParameter(a[i],i);
            }

            cout << "E="<<E<<",  Ep="<<Ep<<",  Em="<<Em<<endl;
            s=h*(Em-Ep)/(2*Ep+2*Em-4*E);

            // TEST
            cout << "s: "<<s<<endl;


            // set new value a
            for(int i=0;i<nParams;i++)
            {
                a[i]+=s*d[i];
                fun->setParameter(a[i],i);
            }


            Newtondiff=fabs(sold-s);
            cout<<"Newtondiff="<<Newtondiff<<endl;
            sold = s;
            Newtoncounter++;
        }




        // calculate -gradient            cout<<"s="<<s<<endl<<"d="<<endl<<d<<endl;

        MC.integrate(fun,H,MC.thermalise(fun,idumadress,parameters),idumadress,parameters);
        here = MC.get_value();
        for(int i=0;i<nParams;i++)
        {
            // move on step forward
            fun->setParameter(a[i]+h,i);
            Rinitial =MC.thermalise(fun,idumadress,parameters);
            MC.integrate(fun,H,Rinitial,idumadress,parameters);
            forward[i]=MC.get_value();

            // set back to initial value
            fun->setParameter(a[i],i);
        }
        gradient = ( -forward + ones(nParams)*here )/h;

        //TEST
        cout<<"Gradient: "<<endl<<ones(nParams)*here<<"-"<<forward<<"/"<<h<<"="<<gradient<<endl;

        resLength = dot(gradient,gradient);
        cout<<"ResLength="<<resLength<<endl;
        if (resLength<resLimit)
        {
            break;
        }

        if (!(counter%nParams))
        {
            t=0.0;
            //TEST
            cout << "outer loop resetet"<<endl;
        }
        else
        {
            t=dot(gradient,gradient-r)/resLength;
            t=max(0.0,t);
            r=gradient;
        }
        d=r+t*d;


        counter++;
    }

    delete fun;
    delete H;
    results.close();

    return 0;
}




// EOF
