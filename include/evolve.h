#ifndef EVOLVE_H
#define EVOLVE_H

#include "objects.h"
#include <stdlib.h>
#include <time.h>
#include <armadillo>


//Parameters regarding evolution

using namespace arma;


class Evolve : public Objects
{
    public:
    int up_or_down,arrive,leave,old_occupation,start,direction,chooser_if_double;
    double gutz_old,gutz_new,eloc;

    void random_hopping();
    double get_K();
    void evolve_A();
    void evolve_iconf();
    int gutz();
    int propose_hop_sistematic();
    int propose_hop_random();
    int metropolis();
    void kels_evolve();
    double calculate_eloc();
};

#endif // EVOLVE_H
