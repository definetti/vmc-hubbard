#include "objects.h"
#include <stdlib.h>
#include <time.h>

/**

A (up and down) are actually the inverse of the matrices made by the occupied rows of a
a is the big matrix wrom which I take the rows
iconf is the system
kelup and keldo are the vectors the places corrisponding to the n-th up or down spin

*/

void Objects::iconf_init_half() //NEEDS nothing PRODUCES iconf. HALF FILLING
{
    int i=0,j=n_electrons/2;
    iconf=zeros<vec>(n_sites);
    do
    {
        i=rand() % n_sites; //piazzo gli spin up (metà)
        if (iconf(i)==0)
        {
            iconf(i)=1;
            j--;
        }
    } while (j > 0);
    for (i=0;i<n_sites;i++) // tutto il resto è spin down
    {
        if(iconf(i)==0)
        {
        iconf(i)=-1;
        }
    }
};

void Objects::iconf_init_gen() //NEEDS nothing PRODUCES iconf. GENERAL
{
    int i=0,j=n_electrons/2;
    iconf=zeros<vec>(n_sites);
    do
    {
        i=rand() % n_sites; //piazzo gli spin up
        if (iconf(i)==0)
        {
            iconf(i)=1;
            j--;
        }
    } while (j > 0);
    j=n_electrons/2;
    do
    {
        i=rand() % n_sites; //piazzo gli spin down
        if (iconf(i)==0)
        {
            iconf(i)=-1;
            j--;
        }
    } while (j > 0);
};

void Objects::kels_init() //NEEDS iconf PRODUCES kelup keldo
{
    int i=0,j=1,h=1;

    kelup=zeros<vec>(n_sites);
    keldo=zeros<vec>(n_sites);

    for(i=0;i<n_sites;i++)
    {
        if(iconf(i)==1 || iconf(i)==2)
        {
            kelup(i)=j;
            j++;
        }
        if(iconf(i)==-1 || iconf(i)==2)
        {
            keldo(i)=h;
            h++;
        }
    }
};

void Objects::a_init() //NEEDS nothing PRODUCES a
{
    int i,j;
    a=zeros<mat>(n_sites,n_sites);
    mat dummy=zeros<mat>(n_sites,n_sites);
    vec eigval; //debug


    for(i=0;i<n_sites;i++)
    {
        for(j=0;j<n_sites;j++)
        {
            if(abs(i-j)==1)
            {
                dummy(i,j)=-t;
            }
        }
    }
    dummy(0,n_sites-1)=-t; //PBC
    dummy(n_sites-1,0)=-t; //PBC
    eig_sym(eigval,a,dummy);
    a.shed_cols(n_electrons/2,n_sites-1);
}

void Objects::A_init() //NEEDS iconf, kelup, keldo,a. PRODUCES A (actually the inverse)
{
    int i;
    Aup=zeros<mat>(n_electrons/2,n_electrons/2);
    Ado=zeros<mat>(n_electrons/2,n_electrons/2);

    for(i=0;i<n_sites;i++)
    {
        if(iconf(i)==1 ||iconf(i)==2)
        {
            Aup.row(kelup(i)-1)=a.row(i);
        }
        if(iconf(i)==-1 || iconf(i)==2)
        {
            Ado.row(keldo(i)-1)=a.row(i);
        }
    }
    Aup=Aup.i();
    Ado=Ado.i();
}

double Objects::kronecker(int a, int b)
{
  if(a==b)
  {
    return 1.0;
  }
  else
  {
    return 0.0;
  }
}
