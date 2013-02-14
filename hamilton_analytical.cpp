#include "hamilton_analytical.h"

hamilton_analytical::hamilton_analytical(Config * parameters) : hamilton(parameters)
{
}


double hamilton_analytical::localEnergy(TrialFct * trialfct, positions * R)
{
    double alpha = trialfct->get_alpha();
    double beta  = trialfct->get_beta();

    // specified for two particles only
    try
    {
        if (nParticles==2)
        {
            double r1 = R->get_r(0);
            double r2 = R->get_r(1);
            double r12= R->get_rr(0,0);
            mat Pos = R->get_pos();
            double r1dotr2 = dot( Pos.col(0) , Pos.col(1) );

            double factor = 1+beta*R->get_rr(0,0);
            double factorsqr = factor*factor;

            double Ee1 = (alpha - Z)*( 1/R->get_r(0) + 1/R->get_r(1)) +1/R->get_rr(0,0) - alpha*alpha;
            double Ee2 = alpha*(r1+r2)/r12 * (1- r1dotr2/(r1*r2) ) -0.5/factorsqr -2/r12 + 2*beta/factor;
            Ee2 /= 2*factorsqr;

            return Ee1 + Ee2;
        }

        else
        {
            throw 2;
        }
    }
    catch (int e)
    {
            cout << "hamilton::analytical_localEnergy is not specified for nParticles !=2." << endl;
    }

}
