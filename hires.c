#include <stdio.h>
#include <math.h>
#define PI 3.14159265358979
#define ABS(a) ((a) < (0) ? (-a) : (a))
#include "jl.h"



#define DIAG1 0
#define MAX(a,b) ((a) >= (b) ? (a) : (b))





/*-------------------------------------------------------*/
/*-------------------------------------------------------*/
/*-------------  HIRES  ----------------------------------*/
/*-------------------------------------------------------*/
int
hires(double *sqr_spec,  double *el, int nwin, int num_freq, double *ares)
{
	int             i, j, k, kpoint;
	float           a;

	for (j = 0; j < num_freq; j++)
		ares[j] = 0.;

	for (i = 0; i < nwin; i++) {
		k = i * num_freq;
		a = 1. / (el[i] * nwin);
		for (j = 0; j < num_freq; j++) {
			kpoint = j + k;
			ares[j] = ares[j] +
				a * ( sqr_spec[kpoint] );
		}
	}

	for (j = 0; j < num_freq; j++) {
		if(ares[j]>0.0) 
                   ares[j] = sqrt(ares[j]);
                  else printf("sqrt problem in hires pos=%d %f\n", j, ares[j]);
	}

	return 1;
}      



