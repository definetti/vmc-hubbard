#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "pp.h"

/**

Does postprocessing over calculated data. Needs number of steps and vector of results

*/

param PP::jackknife (int PASSI, vec v){

	int i,j,k=0,binlength;
	int BINS=PASSI/100;
	param result={0,0};

	binlength=PASSI/BINS;
	vec binsvector = zeros<vec>(BINS);
	vec jack = zeros<vec>(BINS);

	for(i=0;i<PASSI;i=i+binlength) //calculate the average over each bin and puts in a vector

	{
		for(j=i;j<i+binlength;j++)
		{
			binsvector(k)=binsvector(k)+v(j)/binlength;
		}
		k++;
	}

	for(i=0;i<BINS;i++) //calculate the averages without the i-th bin and puts in a vector
	{
		for(j=0;j<BINS;j++)
		{
			if(j!=i)
			{
				jack(i)=jack(i)+(binsvector(j)/(BINS-1));
            }
		}
	}

	for(i=0;i<PASSI;i++) //calculates the mean of the previous quantities and stores in return variable
	{
		result.first+=v(i)/PASSI;
	}

	for(i=0;i<BINS;i++) //calculate jackknife variance and error and stores in return variable
	{
		result.second=result.second+((double)(BINS-1)/(double)(BINS))*(jack(i)-result.first)*(jack(i)-result.first);
	}
	result.second=sqrt(result.second);

	return result;
}
