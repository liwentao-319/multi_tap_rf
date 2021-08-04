/******************************************************************
*
*  modified from multi_tap_spec.c
*  calculate rrf, RCor
*  
*  add gauss window to Z and R spectra before calculate RF spectra 
*  @author liwentao
*******************************************************************/
#include "jl.h"
#include "sac.h"

void zero_pad(float  output[], int start , int olength);
int  multitap(int n, int nwin, double *el,  float npi, double *tapers, double *tapsum);
void jrealft(float data[], unsigned long n, int isign);
int	write_sac(const char	*name,
		SACHEAD		hd,
		const float	*ar
	);

/*---------------------------------------------------------------------------*/
void time_shift(int M0, float *dataf, int klen){
    /*notice t0 = M0*dt and M0 is even,to make sure dataf doesn't change with ff=1/2/dt */
    int ii,jje,jjo;
    float temp,omega;
    for (ii=1;ii<klen/2;++ii){
        jje=ii*2;
        jjo=ii*2+1;
        omega=2*PI*ii/klen;
        temp=(float)cos((double) omega*M0)*dataf[jje] + 
                (float)sin((double) omega*M0)*dataf[jjo];
        dataf[jjo]=(float)cos((double) omega*M0)*dataf[jjo] - 
                (float)sin((double) omega*M0)*dataf[jje];
        dataf[jje]=temp;
    }

}

float gauss_window(float omiga, float alpha){
    float gauss_win;
    gauss_win=-omiga*omiga/(4.*alpha*alpha);
    gauss_win=(float) exp((double)gauss_win);
    return gauss_win;

}
/*---------------------------------------------------------------------------*/
void mt_get_spec(float *series, int inum, int klength, float *amp)
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

/*------------------------------------------------------------------------------*/

void do_mtap_rf_gauss(float *dataz, float *datar, int npts,
	    int nwin, float npi, float dt, float t0, float *rrf,  \
        float *rcorr, int klen, float *noise, float gauss)
{
/*
input:
    dataz = floating point input time signal series in vertical channel
    datar = floating point input time signal series in radial channel 
    npts  = number of points in dataz, datar and datat
    nwin  = number of taper windows to calculate
    npi   = order of the slepian functions
    dt    = sampling interval (time)
    t0    = shift time for datar and datat
    rrf   = spectra output for radial receiver function
    rcorr = correlation in frequency between dataz and datar 
    klen  = number of frequecies calculated (power of 2)
    noise = floating point input time noise series in vertical channel
    npts_n= number of points in noise 
notice:
    rrf, rcorr  must be length of klen
*/

	int             i, j, k, je, jo, M0;
	int             iwin, kk;

	/*************/
	double         *tapsum;	
    double         *lambda, *tapers;  
    float          *ampz,*ampr, *ampn;
    float          ampzz,amprr;
    float          amp_noise,domi;
    float          *dataz_tap, *datar_tap, *noise_tap;
float          *nnn,*zzz, *rrr, *rcorr2;
    float          omiga,window,df=1/dt/klen;
    float          water=0.001;
	/************/
    int             num_freqs;
    int             len_taps, num_freq_tap,num_freq_fft;
	int             kf;
   /**************/ 

#ifdef SAVEAMP
extern SACHEAD HD;
extern char SAVE[128];
char sacout1[128],sacout2[128],sacout3[128],sacout4[128];
strcpy(sacout1,SAVE);strcat(sacout1,".noise.amp");
strcpy(sacout2,SAVE);strcat(sacout2,".zsig.amp");
strcpy(sacout3,SAVE);strcat(sacout3,".rsig.amp");
strcpy(sacout4,SAVE);strcat(sacout4,".rcorr2.amp");
#endif

/*check discretize to to M0, make sure M0=2*j
*/
   M0=t0/dt;
   M0=2*(M0/2);
/*pading zero to noise*/
  /* noise_=(float *) malloc((size_t) npts*sizeof(float));
   if (npts_n>npts) {
       fprintf(stderr,"noise is longer that signal, abort!!!!\n");
       return;
   }
   for (i=0;i<npts;++i){
       noise_[i]=(i<npts_n ? noise[i] : 0);
   }*/
   
/* lambda = vector of eigenvalues   
   tapsum = sum of each taper, saved for use in adaptive weighting  
   tapers =  matrix of slepian tapers, packed in a 1D double array    
*/
    lambda = (double *) malloc((size_t) nwin * sizeof(double));
    tapsum = (double *) malloc((size_t) nwin * sizeof(double));
	len_taps = npts * nwin;
    tapers = (double *) malloc((size_t) len_taps * sizeof(double));
    num_freq_fft = klen*nwin;
/* get a slepian taper  */

	k = multitap(npts, nwin, lambda,  npi, tapers, tapsum);
	
/* apply the taper in the loop.  do this nwin times  */
    ampz = (float *) malloc( (size_t) num_freq_fft*sizeof(float));
    ampr = (float *) malloc( (size_t) num_freq_fft*sizeof(float));
    ampn = (float *) malloc( (size_t) num_freq_fft*sizeof(float));
    
    for (iwin = 0; iwin < nwin; iwin++) {
		kk = iwin * npts;
        kf = iwin * klen;
        dataz_tap = (float *) malloc((size_t) npts*sizeof(float));
        datar_tap = (float *) malloc((size_t) npts*sizeof(float));
        noise_tap = (float *) malloc((size_t) npts*sizeof(float));
		for (j = 0; j < npts; j++){
			dataz_tap[j] = dataz[j] * tapers[kk + j];
            datar_tap[j] = datar[j] * tapers[kk + j];
            noise_tap[j] = noise[j] * tapers[kk + j];
        }

        mt_get_spec(dataz_tap, npts, klen, ampz+kf);  /* calculate the eigenspectrum */
        mt_get_spec(datar_tap, npts, klen, ampr+kf);
        mt_get_spec(noise_tap, npts, klen, ampn+kf);	
    /*time shift*/
        time_shift(M0, ampr+kf, klen);    
        free(dataz_tap);free(datar_tap);free(noise_tap);  
    }  
    
/*calculate reveiver funtion in frequency domain*/
/*check the rrf, ampzz,amprr*/
#ifdef SAVEAMP
zzz=(float *)malloc((size_t)klen/2*sizeof(float));
rrr=(float *)malloc((size_t)klen/2*sizeof(float));
nnn=(float *)malloc((size_t)klen/2*sizeof(float));
rcorr2=(float *)malloc((size_t)klen/2*sizeof(float));
#endif

    for(i=1;i<klen/2;++i){
        rrf[2*i] = 0;
        rrf[2*i+1] = 0;
        omiga=2.*PI*i*df;
        window=gauss_window(omiga,gauss);
        amp_noise = ampn[2*i]*ampn[2*i]+ampn[2*i+1]*ampn[2*i+1];

        /*amp_noise = (float) sqrt((double) amp_noise) ;*/
       
        ampzz=0;amprr=0;
        for (iwin = 0; iwin < nwin; iwin++){
            je = 2*i+klen*iwin;
            jo = je + 1;
            /*amp**2 of channel Z,R and T*/ 
            ampzz += ampz[jo]*ampz[jo] + ampz[je]*ampz[je];
            amprr += ampr[jo]*ampr[jo] + ampr[je]*ampr[je];
            /*R receiver function*/
            rrf[2*i] += ampr[je]*ampz[je] + ampr[jo]*ampz[jo];
            rrf[2*i+1] += ampr[jo]*ampz[je] - ampr[je]*ampz[jo];  
        }

//rfff[i]=rrf[2*i]*rrf[2*i]+rrf[2*i+1]*rrf[2*i+1];

          domi = (float) (sqrt((double) amprr)*sqrt((double) ampzz));
        rcorr[2*i] = rrf[2*i]/domi;
        rcorr[2*i+1] = rrf[2*i+1]/domi;

#ifdef SAVEAMP
zzz[i]=ampzz;
rrr[i]=amprr;
nnn[i]=amp_noise;
rcorr2[i]=rcorr[2*i]*rcorr[2*i]+rcorr[2*i+1]*rcorr[2*i+1];
#endif        
        // rrf[2*i]*=window;
        // rrf[2*i+1]*=window;
          domi = ampzz + amp_noise;
        // domi*=window;
        domi = (domi<=water) ? water : domi;
        rrf[2*i] = rrf[2*i]/domi*window;
        rrf[2*i+1] = rrf[2*i+1]/domi*window;
    }
        /*for i=0,1 responding to f=0,1/(2dt) ,the imaginary part is zero*/
    for (i=0;i<2;++i){
        rrf[i] = 0;
        omiga=2*PI*df*((i==0)?0:klen/2);
        window=gauss_window(omiga,gauss);
        amp_noise = noise[i]*noise[i];
        /*amp_noise = (float) sqrt((double) amp_noise) ;*/

        ampzz=0;amprr=0;
        for (iwin = 0; iwin < 2; iwin++){
            je = i+klen*iwin;
            /*amp**2 of channel Z,R and T*/ 
            ampzz += ampz[je]*ampz[je];
            amprr += ampr[je]*ampr[je];
            /*R receiver function*/
            rrf[i] += ampr[je]*ampz[je];
        }
          domi = (float) (sqrt((double) amprr)*sqrt((double) ampzz));
        rcorr[i] = rrf[i]/domi;

#ifdef SAVEAMP
if (i==0){
    nnn[i]=amp_noise;
    zzz[i]=ampzz;
    rrr[i]=amprr;
    rcorr2[i]=rcorr[i]*rcorr[i];
}
#endif
    //*rrf[i]*=window;
          domi = ampzz + amp_noise;
    //    domi*=window;
        domi = (domi<=water) ? water : domi;
        rrf[i] = rrf[i]/domi*window;
        
    }
#ifdef SAVEAMP
HD.npts=klen/2;
write_sac(sacout1,HD,nnn);
write_sac(sacout2,HD,zzz);
write_sac(sacout3,HD,rrr);
write_sac(sacout4,HD,rcorr2);
#endif
//savedata("zzz.txt",zzz,(unsigned long) klen/2);
//savedata("rrr.txt",rrr,(unsigned long) klen/2);
//savedata("rfff.txt",rfff,(unsigned long) klen/2);

/*  free up memory and return  */
    free(lambda);
    free(tapers);
    free(tapsum);
    free(ampr);
    free(ampz);
    free(ampn);

#ifdef SAVEAMP
    free(zzz);
    free(rrr);
    free(nnn);
    free(rcorr2);
#endif
}
