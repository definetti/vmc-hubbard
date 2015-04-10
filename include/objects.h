#ifndef OBJECTS_H
#define OBJECTS_H

#include <stdlib.h>
#include <armadillo>
#include <time.h>


//Parameters of the system

using namespace arma;


    extern int n_sites;
    extern int n_electrons;
    extern double t;
    extern double U;
    extern double g;

class Objects
{
    public:

    vec kelup;
    vec keldo;
    vec iconf;
    mat a;
    mat Aup;
    mat Ado;

    double kronecker(int a, int b);
    void iconf_init_half();
    void iconf_init_gen();
    void kels_init();
    void a_init();
    void A_init();
};

#endif // OBJECTS_H
