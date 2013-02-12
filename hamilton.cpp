#include "hamilton.h"

// constructor
hamilton::hamilton(Config * parameters)
{
    ndim = parameters->lookup("ndim");
    nParticles = parameters->lookup("nParticles");
    Z = parameters->lookup("Z");
}

// calculate the local energy, 1/phi H phi.
// The Hamiltonian is defined within this function
// \hat{H} = \sum_n (\sum_i \frac12 \del_i^2) -\frac{Z}{r_n} +\sum_{j<i} \frac{1}{r_{ij}}

double hamilton::localEnergy(TrialFct * trialfct, positions * R)
{


    double dens = 0;
    double value = trialfct->getValue(R);

    for (int i=0;i<nParticles;i++)
    {
        // single particle energy
        // potential energy
        dens += -Z/R->get_r(i);
        // kinetic energy
        dens += -0.5*trialfct->getDivGrad(i,R)/value;

        // two paricle part
        for (int j=0; j<i;j++)
        {
            dens += 1/R->get_rr(i-1,j);
        }
    }

    return dens;
}


double hamilton::analytical_localEnergy(TrialFct * trialfct, positions * R)
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
