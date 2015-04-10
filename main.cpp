#include <iostream>
#include <iomanip>      // std::setprecision
#include <fstream>
#include <stdlib.h>
#include <armadillo>
#include <time.h>
#include "objects.h"
#include "evolve.h"
#include "pp.h"


/**

THE PROGRAM SHOULD FUNCTION LIKE THIS:

1. initialize iconf
2. initialize a
3. initialize kels
4. intialize A
5. do metropolis
        5.1. propose hop
        5.2  get K
        GOOD?
            5.2.1 evolve iconf
            5.2.2 evolve A
            5.1.3 calculate Eloc
                5.1.3.1 propose hops
                5.1.3.2 get Ks
            5.2.4 go to point 5
        BAD?
            5.2.1 restore iconf
            5.2.2 Go to point 5

*/


using namespace std;
using namespace arma;

int n_sites;
int n_electrons;
double t;
double U;
double g;


int main(int argc, char *argv[])
{
    ///Variables

    int n_iter,acceptable,counter,ok;
    double controlup, controldo;

    Evolve suzzu; //class evolve
    PP postprocessing; //class pp
    param minima;

    srand (time(NULL));

    counter=0;
    controlup=0;
    controldo=0;
    ok=0;

    ///Parameters initialization

    cerr<< "enter iterations (10^n, n>=3)"<<endl; //>1000 per come ho definito il numero di bin, check pp.cpp
    cin>> n_iter;
    cerr<< "enter number of sites"<<endl;
    cin>> n_sites;
    cerr<<"enter number of electrons"<<endl;
    cin>> n_electrons;


    t=1.0;
    U =2.0;
    g=1.0;

    minima.first=1000; //initialize the variable that stores the minimum of energy, it's obviously enormuos so first step will always be stores and serve as benchmark
    minima.second=1000;
    minima.third=1000;

    ///Vectors initialization

    vec energyvector=zeros<vec>(n_iter);
    do
    {
        if(n_sites!=n_electrons)
        {
            suzzu.iconf_init_gen();
        }
        else
        {
            suzzu.iconf_init_half();
        }
        suzzu.kels_init();
        suzzu.a_init();
        suzzu.A_init();
        controlup=abs(det(suzzu.Aup));
        controldo=abs(det(suzzu.Ado));
    } while (controlup<0.0001 && controldo<0.0001);

    cerr<<"--------------------------"<<endl;
    cerr.precision(0);
    suzzu.iconf.t().raw_print(cerr,"Initial vector");
    cerr.precision(6);
    cerr<<"--------------------------"<<endl;


    ///Cycle for calculation of Eloc
    for(g=0;g<1;g=g+0.1)
    {
        energyvector(0)=suzzu.calculate_eloc();
        //cout <<0<<"   "<<energyvector(0)<< endl;
        counter=1;

        do
        {
            acceptable=suzzu.metropolis(); //accept or reject the move, in case evolves the A
            if(acceptable==1)
            {
                energyvector(counter)=suzzu.calculate_eloc();
                //cout <<counter<<"   "<<energyvector(counter)<< endl;
                counter++;
                ok++;
            }
            else
            {
                energyvector(counter)=energyvector(counter-1);
                //cout <<counter<<"   "<<energyvector(counter)<< endl;
                counter++;
            }
        }while(counter<n_iter);

        postprocessing.energy=postprocessing.jackknife(n_iter,energyvector); //do jackknife

        if(postprocessing.energy.first<minima.first) //store energy minimum
        {
            minima.first=postprocessing.energy.first;
            minima.second=postprocessing.energy.second;
            minima.third=g;
        }

    ///Results

        cout<<g<<"  "<<postprocessing.energy.first<<" "<<postprocessing.energy.second<<endl;

        cerr<<"Done "<<counter<<" iterations for g = "<<g<<endl;
        cerr <<"Acceptance ratio is "<<((double)ok/counter)*100<<" %"<<endl;
        cerr<<""<<endl;
        ok=0;
        counter=0;

    } //end of for

    ///Global results

    cerr<<""<<endl;
    cerr<<"--------------------------"<<endl;
    cerr<<"####### RESULTS #######"<<endl;
    cerr<<"Energy minimum is "<<minima.first<<" ± "<<minima.second<<" for g = "<<minima.third<<endl;
    cerr<<"--------------------------"<<endl;
    return 0;
}
