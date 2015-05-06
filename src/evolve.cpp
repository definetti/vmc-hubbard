#include <stdlib.h>
#include <time.h>
#include "evolve.h"

/**

a is the big matrix from which I take the rows
A (up and down) are actually the inverse of the matrices made by the occupied rows of a
start is the starting site, corresponding to a row of a
leave is the row of A that I leave
arrive is the row of a that I get
up_or_down specifies which spin I am hopping, it's important to store it
gutz-new and gutz_old are the number of doubly-occupied sites befor & after evolution (ok... bad name)
iconf is the system
kelup and keldo are the vectors the places corrisponding to the n-th up or down spin

*/

int Evolve::gutz() //counts number of doubly-occupied sites
{
    int j,n;
    n=0;
    for(j=0;j<n_sites;j++)
    {
        if(iconf(j)==2)
        {
            n++;
        }
    }
    return n;
}

void Evolve::random_hopping() //called from main. calls propose_hop until it proposes a valid hop
{
    int hop_is_possible=0;

    while(hop_is_possible==0)
    {
    	start= rand() % n_sites; //pick the start site for the hopping
        direction=rand() % 2; //pick the hopping direction
        chooser_if_double= rand() % 2; //pick the electron to hop in case site is doubly occupied
    	hop_is_possible=propose_hop_random(); //do the hopping
    }

}

int Evolve::propose_hop_sistematic() //iconf is untouched, produces leave, arrive, gutz_old and gutz_new.
{
    gutz_old=gutz(); //current number of doubly-occupied sites
    gutz_new=gutz_old;
    old_occupation=iconf(start); //store the start site value; il prossimo "if" incasina iconf(start) quindi alla fine devo rimetterlo a posto

    if (iconf(start)==2) //if site is doubly occupied fictiously remove a spin
    {
        gutz_new=gutz_new-1;        //site will not be doubly occupied anymore

        if (chooser_if_double==0) //choose to hop the down
        {
            iconf(start)=-1;
        }
        else      //choose to hop the up
        {
            iconf(start)=1;
        }
    }

    if(direction==0) //left hopping, explicitate arrival sites with PBC
    {
        arrive=(start+n_sites-1) % n_sites;
    }
    else             //right hopping, explicitate arrival sites with PBC
    {
        arrive=(start+n_sites+1) % n_sites;
    }

    up_or_down=iconf(start); //store what I'm hopping (needed in case hop is accepted)

    if((iconf(start)!=iconf(arrive)) && (iconf(arrive)!=2)) //say if the hopping is possible
    {
        if((iconf(start)==1) && (iconf(arrive)==0)) //start site is up, arrive is empty
        {
            leave=kelup(start)-1; //store the value of the column I have left
            iconf(start)=old_occupation;
            return 1;
        }
        if((iconf(start)==1) && (iconf(arrive)==-1)) //start site is up, arrive is down
        {
            leave=kelup(start)-1; //store the value of the column I have left
            gutz_new=gutz_new+1; //there is one more doubly occupied site!
            iconf(start)=old_occupation;
            return 1;
        }
        if((iconf(start)==-1) && (iconf(arrive)==0)) //start site is down, arrive is empty
        {
            leave=keldo(start)-1; //store the value of the column I have left
            iconf(start)=old_occupation;
            return 1;
        }
        if((iconf(start)==-1) && (iconf(arrive)==1)) //start site is down, arrive is up
        {
            leave=keldo(start)-1; //store the value of the column I have left
            gutz_new=gutz_new+1; //there is one more doubly occupied site!
            iconf(start)=old_occupation;
            return 1;
        }
   }
   iconf(start)=old_occupation; //put back the old situation in start site. At the end of this mess iconf has not changed!
   return 0;
};

int Evolve::propose_hop_random() //iconf is untouched, produces leave, arrive, gutz_old and gutz_new.
{
    gutz_old=gutz(); //current number of doubly-occupied sites
    gutz_new=gutz_old;
    old_occupation=iconf(start); //store the start site value; il prossimo "if" incasina iconf(start) quindi alla fine devo rimetterlo a posto

    if(iconf(start)==1 && chooser_if_double==0)
    {
    return 0; //this and the following are crucial. not posing them means probability is not symmetric (singly occupied start site means always hop, doubly means 50% hop!!)
    }

    if(iconf(start)==-1 && chooser_if_double==1)
    {
    return 0;
    }

    if (iconf(start)==2) //if site is doubly occupied fictiously remove a spin
    {
        gutz_new=gutz_new-1;        //site will not be doubly occupied anymore

        if (chooser_if_double==0) //choose to hop the down
        {
            iconf(start)=-1;
        }
        else      //choose to hop the up
        {
            iconf(start)=1;
        }
    }

    if(direction==0) //left hopping, explicitate arrival sites with PBC
    {
        arrive=(start+n_sites-1) % n_sites;
    }
    else             //right hopping, explicitate arrival sites with PBC
    {
        arrive=(start+n_sites+1) % n_sites;
    }

    up_or_down=iconf(start); //store what I'm hopping (needed in case hop is accepted)

    if((iconf(start)!=iconf(arrive)) && (iconf(arrive)!=2)) //say if the hopping is possible
    {
        if((iconf(start)==1) && (iconf(arrive)==0)) //start site is up, arrive is empty
        {
            leave=kelup(start)-1; //store the value of the column I have left
            iconf(start)=old_occupation;
            return 1;
        }
        if((iconf(start)==1) && (iconf(arrive)==-1)) //start site is up, arrive is down
        {
            leave=kelup(start)-1; //store the value of the column I have left
            gutz_new=gutz_new+1; //there is one more doubly occupied site!
            iconf(start)=old_occupation;
            return 1;
        }
        if((iconf(start)==-1) && (iconf(arrive)==0)) //start site is down, arrive is empty
        {
            leave=keldo(start)-1; //store the value of the column I have left
            iconf(start)=old_occupation;
            return 1;
        }
        if((iconf(start)==-1) && (iconf(arrive)==1)) //start site is down, arrive is up
        {
            leave=keldo(start)-1; //store the value of the column I have left
            gutz_new=gutz_new+1; //there is one more doubly occupied site!
            iconf(start)=old_occupation;
            return 1;
        }
   }
   iconf(start)=old_occupation; //put back the old situation in start site. At the end of this mess iconf has not changed!
   return 0;
};

void Evolve::evolve_iconf() //if metropolis is happy evolves iconf using info from propose_hop
{
    int doubleocc=0; //default: the start site is not doubly occupied

    if (iconf(start)==2) //if start site is doubly occupied fictiously remove a spin
    {
        if (up_or_down==1) //choose to hop the down
        {
            iconf(start)=1;
            doubleocc=-1; //remember there was also an up
        }
        else      //choose to hop the up
        {
            iconf(start)=-1;
            doubleocc=1; //remember there was also a down
        }
    }
    if((iconf(start)==1) && (iconf(arrive)==0)) //start site is up, arrive is empty
    {
        iconf(start)=doubleocc; //reput in start site the spins I fictiously took away, if it was doubly occupied
        iconf(arrive)=1; //set the spin in arrival site
        kels_evolve();
        return;
    }
    if((iconf(start)==1) && (iconf(arrive)==-1)) //start site is up, arrive is down
    {
        iconf(start)=doubleocc; //reput in start site the spins I fictiously took away, if it was doubly occupied
        iconf(arrive)=2; //set the spin in arrival site
        kels_evolve();
        return;
    }
   if((iconf(start)==-1) && (iconf(arrive)==0)) //start site is down, arrive is empty
    {
        iconf(start)=doubleocc; //reput in start site the spins I fictiously took away, if it was doubly occupied
        iconf(arrive)=-1; //set the spin in arrival site
        kels_evolve();
        return;
    }
    if((iconf(start)==-1) && (iconf(arrive)==1)) //start site is down, arrive is up
    {
        iconf(start)=doubleocc; //reput in start site the spins I fictiously took away, if it was doubly occupied
        iconf(arrive)=2; //set the spin in arrival site
        kels_evolve();
        return;
    }
};

void Evolve::kels_evolve() //evolvo i kels: sposto l'etichetta dell'elettrone che ha hoppato nel nuovo sito.
{
    if(up_or_down==1)
    {
        kelup(arrive)=kelup(start);
        kelup(start)=0;
    }
    else
    {
        keldo(arrive)=keldo(start);
        keldo(start)=0;
    }
}

double Evolve::get_K() //produces the K factor with info from propose_hop
{
    double K;

    if(up_or_down==1)
    {
        K=dot(a.row(arrive),Aup.col(leave))*exp(-g*(gutz_new-gutz_old));
    }
    else
    {
        K=dot(a.row(arrive),Ado.col(leave))*exp(-g*(gutz_new-gutz_old));
    }
    return K;
};

void Evolve::evolve_A() //evolves A according to info from proposed_hop
{
    int i,j;
    double h=0,r=0;
    mat evolved(n_electrons/2,n_electrons/2);

    if(up_or_down==1)                                       //if I evolve the up matrix
    {
        h=(-1.0)/dot(a.row(arrive),Aup.col(leave));
        for(j=0;j<n_electrons/2;j++)                            //evolve element by element
        {
            r=dot(a.row(arrive),Aup.col(j));
            for(i=0;i<n_electrons/2;i++)
                {
                    evolved(i,j)=Aup(i,j)+h*Aup(i,leave)*(r-kronecker(leave,j));
                }
        }
        Aup=evolved;
    }

    else                                                   //if I evolve the down matrix
    {
        h=(-1.0)/dot(a.row(arrive),Ado.col(leave));
        for(j=0;j<n_electrons/2;j++)                                //evolve element by element
        {
            r=dot(a.row(arrive),Ado.col(j));
            for(i=0;i<n_electrons/2;i++)
                {
                    evolved(i,j)=Ado(i,j)+h*Ado(i,leave)*(r-kronecker(leave,j));
                }
        }
        Ado=evolved;
    }
}

int Evolve::metropolis() //invokes random hopping, gets K and decides if it's ok. If it is evolves iconf and A.
{
    double K,z;

    z= ((double) rand() / (RAND_MAX));
    random_hopping();
    K=get_K();

	if(z < K*K)
	{
        evolve_iconf();
        evolve_A();
        return 1;
    }
    return 0;
};

double Evolve::calculate_eloc() //called from the main if metropolis says ok. Calles propose_hop to get the K. Iconf is untouched.
{
    int response;
    double K,gutz_save;

    gutz_save=gutz_new;//salviamo  il parametro attuale, perchè propose_hop lo cambia
    eloc=0;

    for(start=0;start<n_sites;start++)
    {
        for(direction=0;direction<2;direction++)
        {
            if(iconf(start)==2) //se lo stato è doppiamente occupato devo fare due volte
            {
                for(chooser_if_double=0;chooser_if_double<2;chooser_if_double++)
                {
                    response=propose_hop_sistematic();//cerco i vicini accessibili. Questo dà 0 oppure 1 e modifica iconf (tanto l'abbiamo salvato)
                    if(response==1)
                    {
                        K=get_K(); //come se stessi evolvendo. Serve per eloc
                        eloc=eloc-t*K;
                    }
                }
            }
            else //faccio una volta sola
            {
                response=propose_hop_sistematic();//cerco i vicini accessibili. Questo dà 0 oppure 1 e modifica iconf (tanto l'abbiamo salvato)
                if(response==1)
                {
                    K=get_K(); //come se stessi evolvendo. Serve per eloc
                    eloc=eloc-t*K;
                }
            }
        }
    }
    eloc=eloc+U*gutz_save;//la parte diagonale dell' energia locale dipende dal numero di siti doppiamente occupati.
    return eloc;
}
