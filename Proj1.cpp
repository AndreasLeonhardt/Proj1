#include <iostream>
#include <armadillo>
#include <libconfig.h++>
#include <string>

#include "positions.h"
#include "function.h"
#include "trialfct.h"
#include "trialfct_molecule.h"
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
        fun = new TrialFct_molecule(parameters);
    }

    // create instance of mcInt
    mcInt MC = mcInt(parameters);



    //loop over different parameters alpha, beta and R_0
    int nParams = parameters->lookup("nParameters");
    int parameterIterations = parameters->lookup("parameterIterations");

    vec a=zeros(nParams);


        cout << "..... done."<<endl;


        // PARAMETER OPTIMIZATION ---------------------------------------------------------------------------

        double R_Start = parameters->lookup("R_Start");
        double R0=R_Start;
        double R_Step = parameters->lookup("R_Step");
        int rmax = parameters->lookup("rmax");

        // create names for different files
        char * const samplefile = new char[50];
        char *const pos_samplefile = new char[50];
        const char * samplefilebody = parameters->lookup(("SampleFileName"));
        const char * parameterbody =parameters->lookup("ParameterFile");
        char * const parameterfile = new char [50];
        char * const samplefilePR =new char[50];
        char * const pos_samplefilePR =new char [50];


        for (int r=0;r<rmax;r++)
        {

            cout<<"parameter optimization "<<flush;

            fun->set_R0(R0);

            sprintf(parameterfile,"%s_R%u.txt",parameterbody,r);
            ofstream params;
            params.open(parameterfile);

            a(0)=parameters->lookup("Parameters.[0]");
            a(1)=parameters->lookup("Parameters.[1]");
            fun->setParameter(a);


            // without adaptive stepsize, using 1/i as factor
            for(int i=1;i<parameterIterations+1;i++)
            {
                a -= MC.StatGrad(fun,H,idumadress,nParams,parameters)*(double)(40.0/(40.0+i));
                    fun->setParameter(a);
                    params << a(0)<<"\t"<<a(1)<<endl;
            }
            params.close();

            cout << "..... done."<<endl;


            // INTEGRATION ---------------------------------------------------------------------------------------
            cout <<"integration "<<flush;

            // perform Monte Carlo integration
            // allow for independent Monte Carlo calculations
            // but has to be modified then. pnumber could be a parameter the main() is called with
            // then the blocking and so on should be called speratly.

            sprintf(samplefilePR,"%s%u",samplefilebody,r);
            sprintf(pos_samplefilePR,"%s%u_positions",samplefilebody,r);

            // this loop allows for several random walkers, producing several files with samples.
            // they can be collected. However, the program would have to be changes slightley,
            // for example being called with a different number for 'pnumber'.
            // the parameter 'threats' is more for testing purpose.
            int threats = parameters->lookup("threats");
            for (int pnumber =0;pnumber<threats;pnumber++)
            {
                sprintf(samplefile,"%s_%u.dat",samplefilePR, pnumber);
                sprintf(pos_samplefile,"%s_%u.dat",pos_samplefilePR, pnumber);
                // perform the calculation writing to samplefile
                MC.integrate(fun,H,idumadress,parameters,samplefile,pos_samplefile);
            }
            cout<<"..... done"<<endl;
            // BLOCKING AND WRITING OF RESULTS --------------------------------------------------------------------
            // analyse the result via blocking
            cout<<"blocking "<<flush;

            mat blockingResult = MC.blocking(parameters,samplefilePR,threats).t();

            cout << "...done."<<endl;


            // write results
            cout<<"write results "<<flush;

            ofstream results;
            char * outputfile = new char[50];
            const char * outputfilebody = parameters->lookup("outputfile");
            sprintf(outputfile,"%s_R%u.txt",outputfilebody,r);
            results.open(outputfile);
            results << "Integration points: " << (int) parameters->lookup("nSamples")*(int) parameters->lookup("threats")
                    << "  analytical: " << (int) parameters->lookup("analytical_energy_density") <<endl;
            results << "alpha\tbeta\tR_0\tE\tacceptance_ratio" << endl ;



            results << a(0)<<"\t"<<a(1)<<"\t"<<fun->get_R0()<<"\t"<<MC.get_value()<<"\t"<<MC.get_acceptanceRatio()<<endl<<endl;

            results<<"blocking result:"<<endl
                <<"blocksize\tvariance"<<endl
                <<blockingResult<<endl;

            results.close();

            R0 +=R_Step;

            cout <<"...done"<<endl;

        }



    delete fun;
    delete H;


    // might call a (python)script from here that plots the data in the result file
    // using matplotlib for python for example
    // call the script using the system() function

    char* command=new char[50];
    const char* outputfile = parameters->lookup("outputfile");

    sprintf(command,"python ../Proj1/plot.py %s_R %u %f %f &", outputfile,rmax,R_Start,R_Step);
    cout <<'\n'<<command<<endl;
    system(command);


    return 0;
}




// EOF
