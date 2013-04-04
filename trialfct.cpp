#include "trialfct.h"

// default constructor
TrialFct::TrialFct()
{
    funcParameters[0] = 1.0;
    funcParameters[1] = 1.0;
    stepwidth = 0.001;
    stepwidthsqr = stepwidth*stepwidth;
}

// constructor taking parameters from the config file
TrialFct::TrialFct(Config * parameters) : function(parameters)
{
    funcParameters[0] = parameters->lookup("a_min");
    funcParameters[1] = parameters->lookup("b_min");
    stepwidth = parameters->lookup("stepwidth");
    stepwidthsqr = stepwidth*stepwidth;
}


// calculate the value of the trial function the position given in the variable position.
// Uses hydrogen-like wave functions.

double TrialFct::getValue(positions * R)
{
    double result = 0.0;

    // create Slater matrix
    mat Slaterup  = zeros(nParticleshalf,nParticleshalf);
    mat Slaterdown  = zeros(nParticleshalf,nParticleshalf);



    for (int i=0;i<nParticleshalf;i++)
    {
        for (int j=0;j<nParticleshalf;j++)
        {
            Slaterup(i,j)   = hydrogen(i,j,R);
            Slaterdown(i,j) = hydrogen(i+nParticleshalf,j,R);
        }
    }
    result = det(Slaterup)*det(Slaterdown);


    // jastrow factor
    double jastrow = 0.0;
    double dist;

    for (int i=0;i<nParticleshalf;i++)
    {   // equal spin, factor 1/2
        for (int j=0;j<i;j++)
        {
            dist = R->get_rr(i-1,j);
            jastrow += dist/(2+2*funcParameters[1]*dist);

            dist = R->get_rr(nParticleshalf+i-1,nParticleshalf+j);
            jastrow += dist/(2+2*funcParameters[1]*dist);
        }

        // opposite spin, factor 1/4
        for (int j=nParticleshalf;j<nParticles;j++)
        {
            dist = R->get_rr(j-1,i);
            jastrow += dist/(4+4*funcParameters[1]*dist);
        }

    }

    result *= exp(jastrow);

    return result;
}



// calculate the sum of the numerical second derivatives acting on the trail function ( \nabla^2_i f(x_1,...x_i,...x_n) ),
// derivatives act on the postion of particle n according to the given argument
// uses the position given in the privat variable position.
// using the numerical derivative: f''(x) = 1/h^2 * (f(x+h) + f(x-h) -2*f(x))
double TrialFct::getDivGradOverFct(int particleNumber, positions * R)
{

    double currentValue = getValue(R);
    double value =-2*ndim*currentValue;
    for (int i = 0; i<ndim; i++)
    {
        // move forward: r+h*e_i
        R->step(stepwidth,i,particleNumber);
        value += getValue(R);

        // move backwards: r-h*e_i
        R->step(-2*stepwidth,i,particleNumber);
        value += getValue(R);

        // move to middle
        R->step(stepwidth,i,particleNumber);

    }
    value /=stepwidthsqr*currentValue;
    return value;

}


// calculates numerically the quantum force, defined by 2/(f) * grad(f), where grad(f) is the gradient of f.
// particle number refers to the particle, on whichs position the derivatives act.
vec TrialFct::quantumForce(int particleNumber, positions *R)
{
    vec gradient = zeros(ndim);
    for (int i =0; i<ndim; i++)
    {
        R->step(stepwidth,i,particleNumber);
        gradient(i)+=getValue(R);

        R->step(-2*stepwidth,i,particleNumber);
        gradient(i)-=getValue(R);

        R->step(stepwidth,i,particleNumber);
    }

    gradient/=stepwidth*getValue(R);

    return gradient;

}

//=============================================================================================
// calculate the ratio brut-force
double TrialFct::SlaterRatio(int particleNumber ,positions * Rold,positions * Rnew)
{
    mat Slaternew  = zeros(nParticleshalf,nParticleshalf);
    mat Slaterold  = zeros(nParticleshalf,nParticleshalf);

    if (particleNumber<nParticleshalf)
    {
         for (int i=0;i<nParticleshalf;i++)
         {
             for (int j=0;j<nParticleshalf;j++)
             {
                 Slaternew(i,j) = hydrogen(i,j,Rnew);
                 Slaterold(i,j) = hydrogen(i,j,Rold);
             }
         }
    }
    else
    {
        for (int i=0;i<nParticleshalf;i++)
        {
            for (int j=0;j<nParticleshalf;j++)
            {
                Slaternew(i,j) = hydrogen(i+nParticleshalf,j,Rnew);
                Slaterold(i,j) = hydrogen(i+nParticleshalf,j,Rold);
            }
        }
    }

   return det(Slaternew)/det(Slaterold);

}


void TrialFct::setSlaterinv(positions * R)
{
    // create Slater matrix
    mat Slaterup  = zeros(nParticleshalf,nParticleshalf);
    mat Slaterdown  = zeros(nParticleshalf,nParticleshalf);



    for (int i=0;i<nParticleshalf;i++)
    {
        for (int j=0;j<nParticleshalf;j++)
        {
            Slaterup(i,j)   = hydrogen(i,               j,R);
            Slaterdown(i,j) = hydrogen(i+nParticleshalf,j,R);
        }
    }

    // calculate the inverse of the Slatermatrix for spin up and down
    // respectively. This is used for later updates of the position.
    inverseSlaterDown = inv(Slaterdown);
    inverseSlaterUp   = inv(Slaterup);
}



// Ratio of Jastrow factors, also Brut-Force.
double TrialFct::JastrowRatio(int particleNumber, positions * Rold, positions * Rnew)
{
    // jastrow factor new
    double jastrow_new = 0.0;
    double dist;

    for (int i=0;i<nParticleshalf;i++)
    {   // equal spin, factor 1/2
        for (int j=0;j<i;j++)
        {
            dist = Rnew->get_rr(i-1,j);
            jastrow_new += dist/(2+2*funcParameters[1]*dist);

            dist = Rnew->get_rr(nParticleshalf+i-1,nParticleshalf+j);
            jastrow_new += dist/(2+2*funcParameters[1]*dist);
        }

        // opposite spin, factor 1/4
        for (int j=nParticleshalf;j<nParticles;j++)
        {
            dist = Rnew->get_rr(j-1,i);
            jastrow_new += dist/(4+4*funcParameters[1]*dist);
        }

    }


    // jastrow factor old
    double jastrow_old = 0.0;


    for (int i=0;i<nParticleshalf;i++)
    {   // equal spin, factor 1/2
        for (int j=0;j<i;j++)
        {
            dist = Rold->get_rr(i-1,j);
            jastrow_old += dist/(2+2*funcParameters[1]*dist);

            dist = Rold->get_rr(nParticleshalf+i-1,nParticleshalf+j);
            jastrow_old += dist/(2+2*funcParameters[1]*dist);
        }

        // opposite spin, factor 1/4
        for (int j=nParticleshalf;j<nParticles;j++)
        {
            dist = Rold->get_rr(j-1,i);
            jastrow_old += dist/(4+4*funcParameters[1]*dist);
        }

    }

    return exp(jastrow_new-jastrow_old);
}


void TrialFct::updateSlaterinv(int particleNumber, positions* Rnew, double ratio)
{
    // nothing to do here, since the numerical derivative doesn't use the inverse Slater matrix
}

