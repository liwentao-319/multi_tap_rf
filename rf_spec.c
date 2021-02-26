/***************************************************************************************************
 * multi_tap_rf.c:
 *      Calculate Multi-Taper Receiver functions according to Park and 
 *    Levin, 2000. the program estimating multi-taper spectra is written by
 *    Jonathan M. Lee and Jeffrey Park,referring to the paper in 1994. this
 *    program include subroutines, such as adwait.c, hires.c, ftest.c, jfour1.c,
 *    jrealft.c, jtinvit.c, jtridib.c, mult_tap_spec.c, multitap.c, sigstuff.c. 
 *    those subroutines is included in my program.
 * 
 * Author: 
 *    Wentao lee
 * References:
 *    Lees, J. M. and J. Park. “Multiple-taper spectral analysis: a stand-alone C-subroutine.” 
 *      Computers & Geosciences 21 (1995): 199-236.
 *    Park, J. and V. Levin. “Receiver Functions from Multiple-Taper Spectral Correlation Estimates.” 
 *      Bulletin of the Seismological Society of America 90 (2000): 1507-1520.
 * History:
 *     First Created at 2021/1/20
 ****************************************************************************************************/

#include "jl.h"
#include <assert.h>
#include <string.h>
#include "sac.h"


/*  prototypes  */

int             get_pow_2(int inum); 

void do_mtap_rf(float *dataz, float *datar, int npts,
	    int nwin, float npi, float dt, float t0, float *rrf, \
        float *rcorr,  int klen, float *noise, int  npts_n, float water);

/**********************************************************************************************************/
int main(int argc, char **argv){

    int i,mnwin, error=0;
    int tref, increase;
    int mklen,ok=0;     /*number points for fft*/
    char tracez[128], tracer[128], out1[128], out2[128];
    float  t1, t2, t3, t4, dt, shift, mnpi;
    float *mdataz,*mdatar, *mnoise;
    float *mrrf, *mrcorr;
    float water;
    SACHEAD	hdz,hdr,hds;
    FILE *fout;

	 /*  
     * nwin   ----> the number of summing windows
     * npi    ----> the number of pi-prolate functions, 
     * tref=0 ----> sac header t0 as P arrival that should be write before
     * t1 t2  ----> cut data from t0+t1 to t0+t2 as signal
     * t3 t4  ----> cut data from t0+t3 to t0+t4 as pre-signal noise
     * water  ----> water level to avoid denominator to be zero when deconvoluting
     *  */

	mnpi = 2.5;
	mnwin = 3;
    t1=-5;
    tref=0;
    t2=50;
    t3=-30;
    t4=-5;
    shift=5.;
    water=10e-3;
    /*input argument*/
    ok=0;
    for (i=1; !error && i < argc; i++) {
        if (argv[i][0] == '-') {
            switch(argv[i][1]) {

            case 'S':
                sscanf(&argv[i][2],"%s",tracez);
                ok+=1;
                break;
            case 'R':
                sscanf(&argv[i][2],"%s",tracer);
                ok+=5;
                break;
            case 'W':
                sscanf(&argv[i][2],"%d",&mnwin);
                break;
            case 'P':
                sscanf(&argv[i][2],"%f",&mnpi);
                break;
            case 'T':
                sscanf(&argv[i][2],"%d/%f/%f/%f/%f/%f",&tref,&t1,&t2,&t3,&t4,&shift);
                ok+=7;
                break;
            case 'B':
                sscanf(&argv[i][2],"%f",&water);
                break;
            default:
                error = TRUE;
                break;
        }
        } 
    }

    if (ok!=13||error) {
        fprintf(stderr,"error message:%i\n",ok);
        fprintf(stderr,"usage: %s -Ssrc -Rtrace -W3 -P2 -T0/-5/50/-30/-5/5 \n\
        -S  filename of source time series\n\
        -R  filename of trace time series\n\
        -W  nwin the number of summing windows(optional default 5)\n\
        -P  the number of pi-prolate functions(optional default 3.0)\n\
        -T  format of tref/t1/t2/t3/t4/shift\n\
        -B  water level (optional default 10-4)",argv[0]);
        return -1;
    }
    fprintf(stderr,"tracez:%s\n\
tracer:%s\n\
nwin=%d,npi=%f\n\
tref=%d,t1=%f,t2=%f,t3=%f,t4=%f\n\
water=%f\n\
shift=%f\n",tracez,tracer,mnwin,mnpi,tref,t1,t2,t3,t4,water,shift);
    strcpy(out1,tracer);strcpy(out2,tracer);
    strcat(out1,".rf.spec");
    strcat(out2,".corr.spec");
    /* read sac file*/
    mdataz=read_sac2(tracez,&hdz,tref,t1,t2);
    mdatar=read_sac2(tracer,&hdr,tref,t1,t2);

    if (hdz.npts!=hdr.npts ) {
        fprintf(stderr,"differrent delta bewteen src and trace file\n");
        return -1;
    }
    mnoise=read_sac2(tracez,&hds,tref,t3,t4);
     /*set fft points number */
    
    increase=1;
    mklen = get_pow_2(hdz.npts); 
    mklen=mklen*pow((double) 2,(double) increase);
    printf("channelz delta: %f\n",hdz.delta);
    /* calculate rrf and trf */
    mrrf = (float *) malloc((size_t) mklen*sizeof(float));
    mrcorr = (float *) malloc((size_t) mklen*sizeof(float));
    do_mtap_rf(mdataz, mdatar,hdz.npts, mnwin, mnpi, 
        hdz.delta, shift, mrrf, mrcorr, mklen, mnoise, 
        hds.npts, water);
    /* write mrrf, mtrf, mrcorr and mtcorr*/
    hdz.npts=mklen;
    printf("save rf spectra --> %s\n\
       corr in frequency --> %s\n",out1,out2);
    write_sac(out1,hdz,mrrf);
    write_sac(out2,hdz,mrcorr);
    
    

}