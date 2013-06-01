#ifndef POSITIONS_H
#define POSITIONS_H

#include <armadillo>
#include <libconfig.h++>
#include "lib.h"

using namespace std;
using namespace arma;
using namespace libconfig;

class positions
{
protected:
    int ndim, nParticles;
    mat pos;
    vec r;
    mat rr;

public:
    positions(Config * parameters);
    positions(positions * posit);

    double get_r(int i);
    double get_rr(int i, int j);

    mat get_pos();
    void set_pos(mat NewPositions);
    vec get_singlePos(int particleNumber);
    void set_singlePos(vec NewPosition, int particleNumber);
    void step(double distance, int axis, int Particle);

    int get_ndim();
    int get_nParticles();

};

#endif // POSITIONS_H
