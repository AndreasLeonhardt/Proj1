#include "trialfct.h"

// constructor
TrialFct::TrialFct(Config * parameters, double a, double b )
{
            alpha = a;
            beta = b;
            ndim = parameters->lookup("ndim");
            nParticles = parameters->lookup("nParticles");
            r = vec(nParticles);
            rr = zeros(nParticles,nParticles);


}


// calculate the value of the trial function the position given in the private variable position.
// Here the trial function is defined, currently
// f(r)=e^(-alpha(|r_1|+|r_2|) * e^(|r_1-r_2|/2*(1+beta|r_1-r_2|))
double TrialFct::getValue()
{
    double value = 0;

    for (int i=0;i<nParticles;i++)
    {
        value += -alpha*r[i];
        for (int j=0;j<i;j++)
        {   double dist = rr(i-1,j);
            value += dist/(2+2*beta*dist);
        }
    }
    value = exp( value );

    return value;
}


// overload function for given values
double TrialFct::getValue(mat R,mat RR)
{
    double value = 0;

    for (int i=0;i<nParticles;i++)
    {
        value += -alpha*R[i];
        for (int j=0;j<i;j++)
        {   double dist = RR(i-1,j);
            value += dist/(2+2*beta*dist);
        }
    }
    value = exp( value );

    return value;
}


// calculate the sum of the numerical second derivatives acting on the trail function ( \nabla^2_i f(x_1,...x_i,...x_n) ),
// derivatives act on the postion of particle n according to the given argument
// uses the position given in the privat variable position.
double TrialFct::getDivGrad(int particleNumber)
{
    // stepwidth for numerical differenciation
    double h=0.001;

    double d=0.0;
    double step;

    for (int i = 0; i<ndim; i++)
    {
        // move forward: r+h*e_i
        position(i,particleNumber)+=h;
        updateParticelPosition(position.col(particleNumber),particleNumber);
        step = getValue();


        // move backwards: r-h*e_i
        position(i,particleNumber)-=2*h;
        updateParticelPosition(position.col(particleNumber),particleNumber);
        step += getValue();

        // move to middle
        position(i,particleNumber)+=h;
        updateParticelPosition(position.col(particleNumber),particleNumber);
        d+=step-2*getValue();
    }
    return d/(h*h);
}

// set and get private variables

// change the whole position matrix and update r and rr
void TrialFct::set_position(mat newPosition)
{
    position=newPosition;
    for (int i=0;i<nParticles;i++)
    {
        // length of position vector
        r[i] = norm(position.col(i),2);

        // relative distance
        for (int j=0; j<i;j++)
        {
            rr(i-1,j) = norm(position.col(i)-position.col(j),2);
        }
    }
}

// change a single particle vector and update r and rr
void TrialFct::updateParticelPosition(vec newPosition,int particleNumber)
{
    // write new particle position
    position.col(particleNumber)=newPosition;

    // update r
    r[particleNumber]=norm(position.col(particleNumber),2);

    // update rr
    for (int j=0;j<particleNumber;j++)
    {
        rr(particleNumber-1,j)=norm(position.col(particleNumber)-position.col(j),2);
    }
    for (int i=particleNumber+1; i<nParticles;i++)
    {
        rr(i-1,particleNumber) = norm(position.col(particleNumber)-position.col(i),2);
    }
}


mat TrialFct::get_position()
{
    return position;
}

vec TrialFct::get_r()
{
    return r;
}

mat TrialFct::get_rr()
{
    return rr;
}

void TrialFct::set_alpha(double new_alpha)
{
    alpha = new_alpha;
}

double TrialFct::get_alpha()
{
    return alpha;
}


void TrialFct::set_beta(double new_beta)
{
    beta = new_beta;
}

double TrialFct::get_beta()
{
    return beta;
}
