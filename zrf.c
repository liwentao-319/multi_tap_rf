

#include "jl.h"
#include <string.h>
#include "sac.h"
/*global variable*/
#ifdef SAVEAMP
char SAVE[128];
SACHEAD HD;
#endif
int get_pow_2(int inum);
int strrep(char * str1, char * str2, char * str3);
void do_mtap_rf(float *dataz, float *datar, int npts,
	    int nwin, float npi, float dt, float t0, float *rrf,  \
        float *rcorr, int klen, float *noise);
int mk_dir(char *dir) ;
int	write_sac(const char	*name,
		SACHEAD		hd,
		const float	*ar
	);
void jrealft(float data[], unsigned long n, int isign);

void cos2_oneside_win(float dt, float fc, float *data, int klen);

void mkdirs(char *muldir) ;

void savedata(char *filename,float *data,unsigned long npts);

int main(int argc, char **argv){


    int i,j,ieve,mnwin, error=0,sign,count;
    int tref, increase,isfirst=1;
    int mklen,ok=0,mklen1,imax;     /*number points for fft*/
    float  t1, t2, t3, dt, df, ff, shift, mnpi,fc;
    float *mdataz,*mdatar, *mnoise;
    float *mrrf, *mrcorr, *rf_ave, *averages;
    float cosin,variance,sq_ampc,sq_ampr;
    char saclist[128], zsac[128],sacout[128];
    char datadir[128], outdir[128], paramfile[128],tempz[128],buffer[128];
    SACHEAD	hdz,hdr,hds;
    FILE *fd1;

    #ifdef SAVEAMP
    extern SACHEAD HD;
    extern char SAVE[128];
    #endif

   /*initialize parameters*/
	mnpi = 2.5;
	mnwin = 3;
    tref=0;
    t1=-80;
    t2=-10;
    t3=60;
    shift=5.;
    fc=2.0;
    ok=0;
    strcpy(outdir,"sacout");
    strcpy(datadir,"sacin");
    /*get input parametars from a file */
    for (i=1; !error && i < argc; i++) {
        if (argv[i][0] == '-') {
            switch(argv[i][1]) {

            case 'S':
                sscanf(&argv[i][2],"%s",saclist);
                ok+=1;
                break;
            case 'K':
                sscanf(&argv[i][2],"%d",&mnwin);
                break;
            case 'P':
                sscanf(&argv[i][2],"%f",&mnpi);
                break;
            case 'T':
                sscanf(&argv[i][2],"%d/%f/%f/%f/%f",&tref,&t1,&t2,&t3,&shift);
                ok+=7;
                break;
            case 'F':
                sscanf(&argv[i][2],"%f",&fc);
                ok+=5;
                break;
            case 'I':
                sscanf(&argv[i][2],"%s",datadir);
                break;
            case 'O':
                sscanf(&argv[i][2],"%s",outdir);
                break;
            default:
                error = TRUE;
                break;
        }
        } 
    } 
    if ((t2-t1)!=(t3-t2)){error = TRUE;}  
    /*if error, print the input messages*/ 
    if (ok!=13||error) {
        fprintf(stderr,"error message:%i\n",ok);
        fprintf(stderr,"usage: %s -Ssaclist  -K3 -P2 -T0/-5/50/-30/-5/5 -F2.0  \n\
        -S  filename of source time series\n\
        -K  nwin the number of summing windows(optional default 5)\n\
        -P  the number of pi-prolate functions(optional default 3.0)\n\
        -T  format of tref/t1/t2/t3/shift ( require t2-t1=t3-t2)\n\
        -I  directory to load sac files\n\
        -F  cutoff frequency\n\
        -O  directory to save the results\n",argv[0]);
        return -1;
    }

     /*print the parameters to screen*/    
    fprintf(stderr,"saclistfile:%s\n\
nwin=%d,npi=%f\n\
tref=%d,t1=%f,t2=%f,t3=%f\n\
shift=%f\n\
fc=%f\n\
datadir:%s\n\
output:%s\n",saclist,mnwin,mnpi,tref,t1,t2,t3,shift,fc,datadir,outdir);

    mkdirs(outdir);
   // strcpy(saclist,"testlist.txt");
    fd1=fopen(saclist, "r");
    if (fd1==NULL) {
        fprintf(stderr,"can't open file %s\n",saclist);
        return -1;
    }


    ieve=0;
    while(!feof(fd1)){
        fscanf(fd1,"%s",tempz);
        ieve+=1;
        strcpy(zsac,datadir);
        strcat(zsac,"/");
        strcat(zsac,tempz);
        /*user-self definite*/
        fprintf(stderr, "ievent:%d\n   zfile:%s\n  ",ieve,zsac);
        mdataz=read_sac2(zsac,&hdz,tref,t2,t3);     
        mnoise=read_sac2(zsac,&hds,tref,t1,t2);
        increase=1;
        mklen = get_pow_2(hdz.npts); 
        mklen=mklen*pow((double) 2,(double) increase);
        if (isfirst==1){
            mklen1=mklen;dt=hdz.delta;
            df=1/dt/mklen1;
            imax=fc/df;
            mrrf = (float *) malloc((size_t) mklen1*sizeof(float));
            mrcorr = (float *) malloc((size_t) mklen1*sizeof(float)); 
            rf_ave = (float *) malloc((size_t) mklen1*sizeof(float));
            isfirst=0;
        }
        if (mklen!=mklen1){
            fprintf(stderr,"different points of input sac:%s\n",zsac);
            return -1;
        }
        #ifdef SAVEAMP
        HD=hdz;
        strcpy(SAVE,outdir);
        strcat(SAVE,"/");
        strcat(SAVE,tempz);
        #endif

        /****************************************************/
        do_mtap_rf(mdataz, mdataz,hdz.npts, mnwin, mnpi, 
            hdz.delta, shift, mrrf, mrcorr, mklen, mnoise);
        
 
        /********************************************************************/
        /*invert the Fourier transform for each rf*/ 
        /**************************************************************/        
        cos2_oneside_win(dt,fc,mrrf,mklen);
        sign=-1;
        jrealft(mrrf-1,(unsigned long)mklen1,sign);
        
        /*  mult by nnf/nf to compensate for losing high freqs
        mult by 2 if cosine**2 taper is applied to spectral RF*/
        for (i=0;i<mklen1;++i) mrrf[i]=mrrf[i]*2/dt/fc/mklen1;
        
        /*save each rf to sac file*/
        hdz.npts=mklen1;
        strcpy(sacout,outdir);
        strcat(sacout,"/");
        strcat(sacout,tempz);
        strcat(sacout,".zrf.time");
        write_sac(sacout,hdz,mrrf);  
    }
    return 0;

}
