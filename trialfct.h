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
    double getDivGradOverFct(int particleNumber, positions * R);
    // calculate the quantum force defined by grad(f)/f of the function numerically
    vec quantumForce(int particleNumber, positions *R);

    double SlaterRatio(int particleNumber ,positions * Rold,positions * Rnew);
    // double JastrowRatio(i,Rold,Rnew);//TODO


    void setSlaterinv(positions * R);
    void updateSlaterinv(int particleNumber, positions* Rnew, double ratio);


    double hydrogen(int particleNumber, int orbital, positions * R);
    vec gradhydrogen(int particleNumber, int orbital, positions *R);
    double divgradhydrogen(int particleNumber, int orbital, positions* R);


};

#endif // TRIALFCT_H

