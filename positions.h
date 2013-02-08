#ifndef POSITIONS_H
#define POSITIONS_H

#include <armadillo>
#include <libconfig.h++>

using namespace std;
using namespace arma;
using namespace libconfig;

class positions
{
    int ndim, nParticles;
    mat pos;
    mat r;
    mat rr;

public:
    positions(Config * parameters);

    double get_r(int i);
    double get_rr(int i, int j);

    mat get_pos();
    void set_pos(mat NewPositions);
    void set_singlePos(vec NewPosition, int particleNumber);
    void step(double distance, int axis, int Particle);

};

#endif // POSITIONS_H
