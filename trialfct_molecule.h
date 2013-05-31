#ifndef TRIALFCT_ANALYTICAL_H
#define TRIALFCT_ANALYTICAL_H

#include <armadillo>
#include <libconfig.h++>

#include "positions.h"
#include "function.h"

using namespace std;
using namespace arma;
using namespace libconfig;


class TrialFct_molecule : public function
{


public:
    // constructor, takes the above parameters
    TrialFct_molecule(Config * parameters);
    // default constructor without parameters
    TrialFct_molecule();

    // calculate the value of the trial function at a certain position
    double getValue(positions * R);
    // calculate the sum of the second derivatives acting on the trial function
    double getDivGradOverFct(int particleNumber, positions * R);
    // calculate the quantum force defined by grad(f)/f of the function numerically
    vec quantumForce(int particleNumber, positions *R);


    double SlaterRatio(int particleNumber ,positions * Rold,positions * Rnew);
    double JastrowRatio(int particleNumber, positions * Rold, positions * Rnew);
    double ParamDerivativeOverFct(positions *R, int parameterNumber);

    void setSlaterinv(positions * R);
    void updateSlaterinv(int particleNumber, positions* Rnew, double ratio);

    vec GradJastrow(int particleNumber, positions * R);
    vec GradSlater(int particleNumber, positions * R);

};

#endif // TRIALFCT_ANALYTICAL_H
