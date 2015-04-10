#ifndef PP_H
#define PP_H

#include <stdlib.h>
#include <armadillo>
#include <time.h>

using namespace arma;

struct param                        //a type consisting in three variables, which will be useful
{
	double first;
	double second;
	double third;
};

class PP
{
    public:
    param energy;

    param jackknife (int PASSI, vec v); //calculates mean values and errors
};

#endif // PP_H
