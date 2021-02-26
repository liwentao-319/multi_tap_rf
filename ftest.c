#include "jl.h"


void 
get_F_values(double *sr, double *si, int nf, int nwin, float *Fvalue, double *b)
{
	/*    input:
             sr si are the  real and imaginary parts of the eigenspectra 
              nf = number of frequency points
              nwin = number of taper windows
	     b is fft of slepian eigentapers at zero freq 
             Cee contains line frequency estimates and f-test parameter
             return:
              Fvalue

           see equation 13.10 in Thomson (1982)
           or page 499 in Percival and Walden (1993)
	 */
	double          sum, sumr, sumi, sum2;
	int             i, j, k;
	double         *Ceer, *Ceei;
	sum = 0.;

	Ceer = (double *) malloc((size_t) nf * sizeof(double));
	Ceei = (double *) malloc((size_t) nf * sizeof(double));

	


	for (i = 0; i < nwin; i++) {
                
		sum = sum + b[i] * b[i];
	}
	for (i = 0; i < nf; i++) {
		Ceer[i] = 0.;
		Ceei[i] = 0.;
		for (j = 0; j < nwin; j++) {
			k = i + j * nf;
			Ceer[i] = Ceer[i] + sr[k] * b[j];
			Ceei[i] = Ceei[i] + si[k] * b[j];
		}
		Ceer[i] = Ceer[i] / sum;
		Ceei[i] = Ceei[i] / sum;
		sum2 = 0.;
		for (j = 0; j < nwin; j++) {
			k = i + j * nf;
			sumr = sr[k] - Ceer[i] * b[j];
			sumi = si[k] - Ceei[i] * b[j];
			sum2 = sum2 + sumr * sumr + sumi * sumi;
		}
		Fvalue[i] = (float) (nwin - 1) * (SQR(Ceei[i]) + SQR(Ceer[i])) * sum / sum2;
	}
          free(Ceei);
          free(Ceer);
	return;
}
