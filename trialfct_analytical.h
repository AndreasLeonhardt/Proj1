#ifndef TRIALFCT_ANALYTICAL_H
#define TRIALFCT_ANALYTICAL_H

#include <armadillo>
#include <libconfig.h++>

#include "positions.h"
#include "function.h"

using namespace std;
using namespace arma;
using namespace libconfig;


class TrialFct_analytical : public function
{


    // parameters of the trial function
    double alpha, beta;



public:
    // constructor, takes the above parameters
    TrialFct_analytical(Config * parameters);
    // default constructor without parameters
    TrialFct_analytical();

    // calculate the value of the trial function at a certain position
    double getValue(positions * R);


    // calculate the sum of the second derivatives acting on the trial function
    double getDivGrad(int particleNumber, positions * R);
    // calculate the quantum force defined by grad(f)/f of the function numerically
    vec quantumForce(int particleNumber, positions *R);

    void set_alpha(double new_alpha);
    double get_alpha();

    void set_beta(double new_beta);
    double get_beta();





};

#endif // TRIALFCT_ANALYTICAL_H
