#include "jl.h"


#define perr(x,y)  (fprintf(stderr, x , y))
#define prbl (fprintf(stderr,"\n"))


void zero_pad(float  output[], int start , int olength);
int 
adwait(double *tapers,  double *dcf,
       double *el, int nwin, int nf, double *ares, double *degf, double avar);
int
hires(double *spreal,  double *el, int nwin, int nf, double *ares);

int  multitap(int n, int nwin, double *el,  float npi, double *tapers, double *tapsum);

void  get_F_values(double *sr, double *si, int nf, int nwin,float *Fvalue, double *b);

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

void 
mt_get_spec(float *series, int inum, int klength, float *amp)
{
/*    series = input time series
      inum   = length of time series
      klength = number of elements in power spectrum (a power of 2)
      amp = returned power spectrum
*/

	int             i, j, isign = 1;

	unsigned long   nn;
	float           tsv;


    	nn = klength;

        

	/* copy amp onto series and apply zero padding to  klength */
       
	for (i = 0; i < inum; i++) {
		
		amp[i] = series[i];
	
	}

       zero_pad(amp, inum, klength);


/*  Fast Fourier Transform Routine:  here we are using the Numerical Recipes
     routine jrealft which returns the fft in the 1-D input array
     packed as pairs of real numbers.
     The jrealft routine requires the input array to start at index=1
      so we must decrement the index of amp
*/


 
          jrealft(amp-1, nn, isign);


}

void 
do_mtap_spec(float *data, int npoints, int kind,
	    int nwin, float npi, int inorm, float dt, float *ospec, float *dof, float *Fvalues, int klen)
{
/*
     data = floating point input time series
    npoints = number of points in data
     kind = flag for choosing hires or adaptive weighting coefficients
     nwin = number of taper windows to calculate
     npi = order of the slepian functions
     inorm = flag for choice of normalization
     dt = sampling interval (time)
     ospec = output spctrum
    dof = degrees of freedom at each frequency
    Fvalues = Ftest value at each frequency estimate
    klen = number of frequecies calculated (power of 2)

*/

	int             i, j, k;
	double         *lambda, *tapers;
	long            len, longlen;
	float          *xt;
	FILE           *fopen(), *inf, *tapfile;
        FILE            *dof_file;

	int             logg;
	int             nn;
	float          *b;
	int             iwin, kk;

	/*************/
	double          anrm, norm;
        double            *ReSpec, *ImSpec;
	double         *sqr_spec,  *amu;
	float          *amp, *fv;
	double          avamp, temp, sqramp;
	double          sum, *tapsum;
	/************/
         int num_freqs;
         int len_taps, num_freq_tap;
   
	double         *dcf, *degf, avar;
	int             n1, n2, kf;
	int             flag;
	int             one = 1;

         double tem1, tem2;

/* lambda = vector of eigenvalues   
   tapsum = sum of each taper, saved for use in adaptive weighting  
   tapers =  matrix of slepian tapers, packed in a 1D double array    
*/
lambda = (double *) malloc((size_t) nwin * sizeof(double));
tapsum = (double *) malloc((size_t) nwin * sizeof(double));


        
	len_taps = npoints * nwin;

tapers = (double *) malloc((size_t) len_taps * sizeof(double));

             num_freqs = 1+klen/2;
             num_freq_tap = num_freqs*nwin;



       /* get a slepian taper  */

	k = multitap(npoints, nwin, lambda,  npi, tapers, tapsum);
#if 0
        /* print out tapers for curiosity  */
        for(i=0; i<npoints; i++){
         for(j=0; j<nwin; j++)fprintf(stderr,"%d %15.10f ",i,tapers[i+j*npoints]);
          prbl;
          }
#endif



	
     /* choose normalization based on inorm flag  */

	anrm = 1.;

	switch (inorm) {
	case 1:
		anrm = npoints;
		break;
	case 2:
		anrm = 1 / dt;
		break;
	case 3:
		anrm = sqrt((double) npoints);
		break;
	default:
		anrm = 1.;
		break;
	}

	
	/* apply the taper in the loop.  do this nwin times  */


amu = (double *) malloc((size_t) num_freqs * sizeof(double));
sqr_spec = (double *) malloc((size_t) num_freq_tap * sizeof(double));
ReSpec = (double *) malloc((size_t) num_freq_tap * sizeof(double));
ImSpec = (double *) malloc((size_t) num_freq_tap * sizeof(double));


	


	for (iwin = 0; iwin < nwin; iwin++) {
		kk = iwin * npoints;
                kf = iwin * num_freqs;

b = (float *) malloc((size_t) npoints * sizeof(float));

		for (j = 0; j < npoints; j++)
			b[j] = data[j] * tapers[kk + j];   /*  application of  iwin-th taper   */




amp = (float *) malloc((size_t) klen * sizeof(float));

		
	


 
		mt_get_spec(b, npoints, klen, amp);  /* calculate the eigenspectrum */

	        free(b);
          
		
		sum = 0.0;


/* get spectrum from real fourier transform    */


          norm = 1.0/(anrm*anrm);




            for(i=1; i<num_freqs-1; i++){
       if(2*i+1 > klen) fprintf(stderr,"error in index\n");
       if(i+kf > num_freq_tap ) fprintf(stderr,"error in index\n");

            sqramp = SQR(amp[2*i+1])+SQR(amp[2*i]);

            ReSpec[i+kf] = amp[2*i];
            ImSpec[i+kf] = amp[2*i+1];



            sqr_spec[i+kf] =    norm*(sqramp);

             sum += sqramp;
            }
          sqr_spec[0+kf] = norm*SQR(fabs(amp[0]));
          sqr_spec[num_freqs-1+kf] = norm*SQR(fabs(amp[1]));

            ReSpec[0+kf] = amp[0];
            ImSpec[0+kf] = 0.0;

            ReSpec[num_freqs-1+kf] = amp[1];
            ImSpec[num_freqs-1+kf] = 0.0;

             sum += sqr_spec[0+kf] + sqr_spec[num_freqs-1+kf];

        if(num_freqs-1+kf>num_freq_tap )fprintf(stderr,"error in index\n");

		temp = sum / (double) num_freqs;
		if (temp > 0.0)
			avamp = sqrt(temp) / anrm;
		else {
			avamp = 0.0;
			 /* fprintf(stderr," avamp = 0.0! \n"); */ 
		}


		free(amp);

	}


fv = (float *) malloc((size_t) num_freqs * sizeof(float));

	
        /* choice of hi-res or adaptive weighting for spectra    */

	switch (kind) {
	case 1:
	
		hires(sqr_spec,  lambda, nwin, num_freqs, amu);
	        get_F_values(ReSpec, ImSpec, num_freqs, nwin, fv, tapsum);
 

 
           for (i = 0; i < num_freqs; i++) {
		ospec[i] =amu[i];
                 dof[i] = nwin-1;
                 Fvalues[i] = fv[i];
                  }

		break;

	case 2:


		/* get avar = variance*/

		n1 = 0;
		n2 = npoints;


                  avar = 0.0;

		for (i = n1; i < n2; i++)
			avar += (data[i]) * (data[i]);


		switch (inorm) {
		case 1:
			avar = avar / (npoints * npoints);
			break;

		case 2:
			avar = avar * dt * dt;
			break;

		case 3:

			avar = avar / npoints;
			break;

		default:
			break;
		}

		 
	dcf = (double *) malloc((size_t) num_freq_tap * sizeof(double));
	degf = (double *) malloc((size_t) num_freqs * sizeof(double));


	
	


		adwait(sqr_spec, dcf, lambda, nwin, num_freqs, amu, degf, avar);

                get_F_values(ReSpec, ImSpec, num_freqs, nwin, fv, tapsum);

#if 0
           /* dump out the degrees of freedom to a file for later inspection  */
              	if ((dof_file = fopen("dof_file", "w")) == NULL) {
		fprintf(stderr, "dof unable to open\n");
		return;}
               	for (i = 0; i < num_freqs; i++) {
	                fprintf(dof_file,"%f\n",degf[i]);
                       	}

                   fclose(dof_file);
#endif

                 /* rap up   */

           for (i = 0; i < num_freqs; i++) {
		ospec[i] =amu[i];
                 dof[i] = degf[i];
                 Fvalues[i] = fv[i];
                  }


                
		free(dcf);
		free(degf);
		free(fv);


		break;
	}

/*  free up memory and return  */

        free(amu);

	free(sqr_spec);
	free(ReSpec);

	free(ImSpec);

	free(lambda);

	free(tapers);


}
