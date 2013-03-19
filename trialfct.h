#ifndef TRIALFCT_H
#define TRIALFCT_H

#include <armadillo>
#include <libconfig.h++>

#include "positions.h"
#include "function.h"

using namespace std;
using namespace arma;
using namespace libconfig;


class TrialFct : public function
{
    // parameters of the trial function
    //double parameters[2];

    // parameters for numerical differentiation
   // double stepwidth, stepwidthsqr;




public:
    // constructor, takes the above parameters
    TrialFct(Config * parameters);
    // default constructor without parameters
    TrialFct();

    // calculate the value of the trial function at a certain position
    double getValue(positions * R);


    // calculate the sum of the second derivatives acting on the trial function
    double getDivGrad(int particleNumber, positions * R);
    // calculate the quantum force defined by grad(f)/f of the function numerically
    vec quantumForce(int particleNumber, positions *R);


    void setParameter(double newParameter, int parameterNumber);
    double getParameter(int parameterNumber);
    double hydrogen(int particleNumber, int orbital, positions * R);

    void set_stepwidth(double new_stepwidth);


};

#endif // TRIALFCT_H

