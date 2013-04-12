#ifndef FUNCTION_H
#define FUNCTION_H


#include <armadillo>
#include <libconfig.h++>

#include "positions.h"

using namespace std;
using namespace arma;
using namespace libconfig;


class function
{
protected:

    // more general parameters (position vector, number of spacial dimensions, number of particles.
    int ndim, nParticles, nParticleshalf,nParams;
    // number of parameters for the function needs obviously to be fixed here.
    // in the case of varying parameters setting this to the maximal value would be
    // the quick and dirty solution. So far, 2 parameter seem to be sufficient.
    // else, one could use maybe an armadillo vector
    double funcParameters[2];
    // parameters for numerical derivation
    double stepwidth, stepwidthsqr;

    mat inverseSlaterUp;
    mat inverseSlaterDown;

public:
    function();
    function(Config * parameters);


    // calculate values of the function and its derivatives at the positon &R,
    // derivatives acting on the coordinates of particle particleNumber

    // value of the function
    virtual double getValue(positions * R)=0;

    // sum of the second derivatives, div(grad(f))
    virtual double getDivGradOverFct(int particleNumber, positions * R)=0;

    // quantum force defined by g*grad(f)/f
    // g to be specified (1/2 for fermionic wave functions
    virtual vec quantumForce(int particleNumber, positions *R)=0;


    virtual double SlaterRatio(int particleNumber ,positions * Rold,positions * Rnew)=0;
    virtual double JastrowRatio(int particleNumber, positions * Rold, positions * Rnew)=0;


    virtual void setSlaterinv(positions * R)=0;
    virtual void updateSlaterinv(int particleNumber, positions* Rnew, double ratio)=0;
    virtual double ParamDerivativeOverFct(positions *R,int parameterNumber)=0;

    double hydrogen(int particleNumber, int orbital, positions * R);
    vec gradhydrogen(int particleNumber, int orbital, positions *R);
    double divgradhydrogen(int particleNumber, int orbital, positions* R);
    double dhydrogenda(int particleNumber, int orbital,positions* R);




    // get parameter for numerical derivatives.
    // this would be nicer if just available in the derived class with "_numerical"
    // However, this variables seem to be needed here as well, in order to set them to proper values.
    double get_stepwidth();
    double get_stepwidthsqr();

    void setParameter(double newParameter, int parameterNumber);
    void setParameter(vec newParameters);
    double getParameter(int parameterNumber);

    void setndim(int newvalue);
    int getndim();

    void setnParticles(int numberofParticles);
    int getnParticles();

    mat getinvslatermatrix(int particleNumber);
    void setSlaterinv(int particleNumber, mat newSlaterInv);




};

#endif // FUNCTION_H
