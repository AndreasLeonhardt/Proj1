#include <iostream>
#include<armadillo>

#include "trialfct.h"

using namespace std;
using namespace arma;



int ndim = 3;
int nParticles = 2;



int main()
{   // seed for randu
    srand(41);
    // the position of particle i is stored in the i-th column.
    mat positions =randu(ndim,nParticles);

    // parameters of trial wave function
    double alpha = 0.5;
    double beta = 1;

    // create instance of TrialFct
    TrialFct fun = TrialFct(alpha, beta, ndim, nParticles);
    fun.set_position(positions);

    cout << "f(r1,r2)=" << fun.getValue() <<"  D_^2 f(r1,r2)=" << fun.getDivGrad(0) << endl;

    return 0;
}



// EOF
