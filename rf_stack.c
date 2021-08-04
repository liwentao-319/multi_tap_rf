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
void do_mtap_rf_gauss(float *dataz, float *datar, int npts,
	    int nwin, float npi, float dt, float t0, float *rrf,  
        float *rcorr, int klen, float *noise,float gauss);
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
    int i,j,neve,mnwin,error=0,sign,count,nstack,stack_t=0;
    int tref, increase,isfirst=1,error1=0,*counts_stack;
    int mklen,ok=0,mklen1,imax;     /*number points for fft*/
    float  t1, t2, t3, dt, df, ff, shift, mnpi,gauss,minstack,maxstack,dstack,stack_range,stack;
    float *mdataz,*mdatar, *mnoise;
    float *mrrf, *mrcorr, *rf_ave, *averages;
    float cosin,variance,sq_ampc,sq_ampr;
    char saclist[128], zsac[128], rsac[128],sacout[256],schannel[12],rchannel[12];
    char datadir[128], outdir[128], paramfile[128],tempz[128],tempr[128],buffer[128];
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
    ok=0;
    gauss=3.0;
    strcpy(rchannel,".R.");/*channel to deconvolute. */
    strcpy(schannel,".Z.");/*soure time series used to convolute rchannel*/
    strcpy(outdir,"sacout");
    strcpy(datadir,"sacin");
    minstack=30;
    maxstack=90;
    dstack=10;
    stack_range=5;
    
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
            case 'g':
                sscanf(&argv[i][2],"%f/%f/%f/%f/%d",&minstack,&maxstack,&dstack,&stack_range,&stack_t);
                break;
            case 'G':
                sscanf(&argv[i][2],"%f",&gauss);
                ok+=5;
                break;
            case 'I':
                sscanf(&argv[i][2],"%s",datadir);
                break;
            case 'O':
                sscanf(&argv[i][2],"%s",outdir);
                break;
            case 'R':
                sscanf(&argv[i][2],"%s",rchannel);//channel marker for numerator of deconvolution
                break;
            case 'Z':
                sscanf(&argv[i][2],"%s",schannel);//channel marker for denominator of deconvolution
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
        fprintf(stderr,"usage: %s -Ssaclist  -K3 -P2 -T0/-5/50/-30/-5/5 -G2.0 -Idirin -Odirout \n\
        -S  filename of source time series\n\
        -K  nwin the number of summing windows(optional default 5)\n\
        -P  the number of pi-prolate functions(optional default 3.0)\n\
        -T  format of tref/t1/t2/t3/shift ( require t2-t1=t3-t2)\n\
        -I  directory to load sac files\n\
        -G  gauss parameter\n\
        -O  directory to save the results\n\
        -g  minstack/maxstack/dstack/stack_type 0-->stack|1-->baz\n",argv[0]);

        return -1;
    }

     /*print the parameters to screen*/    
    fprintf(stderr,"saclistfile:%s\n\
nwin=%d,npi=%f\n\
tref=%d,t1=%f,t2=%f,t3=%f\n\
shift=%f\n\
gauss=%f\n\
datadir:%s\n\
output:%s\n\
min:%f   max:%f   increment:%f  sum_range:%f\n",saclist,mnwin,mnpi,tref,t1,t2,t3,shift,gauss,datadir,\
                        outdir,minstack,maxstack,dstack,stack_range);
    nstack=(maxstack-minstack)/dstack+1;
    mkdirs(outdir);
   // strcpy(saclist,"testlist.txt");
    fd1=fopen(saclist, "r");
    if (fd1==NULL) {
        fprintf(stderr,"can't open file %s\n",saclist);
        return -1;
    }

    
    neve=0;
    while(!feof(fd1)){
        fscanf(fd1,"%s",tempz);
        neve+=1;
        strcpy(zsac,datadir);
        strcat(zsac,"/");
        strcat(zsac,tempz);
        strcpy(rsac,zsac);
        if ((count = strrep(rsac,schannel,rchannel))==0){
            fprintf(stderr,"error in replace str:%s\n",rsac);
            return -1;
        };/*user-self definite*/

        strcpy(tempr,tempz);

        if ((count = strrep(tempr,schannel,rchannel))==0){
            fprintf(stderr,"error in replace str:%s\n",tempr);
            return -1;
        };
        /*user-self definite*/
        mdataz=read_sac2(zsac,&hdz,tref,t2,t3,&error1);
        if (error1==-1){
            fprintf(stderr,"Can't open sacfile:%s, for signal\n",zsac);
            error1=0;
            continue;
        } 
        mdatar=read_sac2(rsac,&hdr,tref,t2,t3,&error1);
        if (error1==-1){
            fprintf(stderr,"Can't open sacfile:%s, for signal\n",rsac);
            error1=0;
            continue;
        } 
        if (hdz.npts!=hdr.npts ) {
            fprintf(stderr,"differrent delta bewteen src and trace file\n");
            return -1;
        }      
        mnoise=read_sac2(zsac,&hds,tref,t1,t2,&error1);
        if (error1==-1){
            fprintf(stderr,"Can't open sacfile:%s, for noise\n",zsac);
            error1=0;
            continue;
        } 


        /*****************************************************/
        increase=1;
        mklen = get_pow_2(hdz.npts); 
        mklen=mklen*pow((double) 2,(double) increase);
        /*for the first loop allocate memory*/
        if (isfirst==1){
            mklen1=mklen;dt=hdz.delta;
            df=1/dt/mklen1;
            mrrf = (float *) malloc((size_t) mklen1*sizeof(float));
            mrcorr = (float *) malloc((size_t) mklen1*sizeof(float)); 
            rf_ave = (float *) calloc( nstack,mklen1*sizeof(float));
            averages = (float *) calloc( nstack,mklen1*sizeof(float));
            counts_stack = (int *)malloc( nstack*sizeof(int));
            for (i=0;i<nstack;++i) counts_stack[i]=0;
            isfirst=0;
        }
        /*********************************************************/
        if (mklen!=mklen1){
            fprintf(stderr,"different points of input sac:%s\n",zsac);
            return -1;
        }
        #ifdef SAVEAMP
        HD=hdz;
        strcpy(SAVE,outdir);
        strcat(SAVE,"/");
        strcat(SAVE,tempr);
        #endif

        /*calculate spectra of both source time series and input signal time series by multi-taper*/
        do_mtap_rf_gauss(mdataz, mdatar,hdz.npts, mnwin, mnpi, 
            hdz.delta, shift, mrrf, mrcorr, mklen, mnoise,gauss);
        
        /**************************************************************/ 

        /*chose stack type according to gcarc(0) or baz(1)*/
        if (stack_t==0){
            stack=hdz.gcarc;
        }
        else if (stack_t==1){

            stack=hdz.baz;
        }
        else{
            fprintf(stderr,"invalid input for stack_type");
            return -1;
        }
        if (stack<minstack || stack>maxstack) continue;
        /*loop though all range of stack points*/
        fprintf(stderr,"eventid:%d---%s\n  stack at:",neve,tempz);
        for (j=0;j<nstack;++j){

            if ((stack<minstack+j*dstack-stack_range) || (stack>minstack+j*dstack+stack_range) ) continue;
            fprintf(stderr,"%f   ", minstack+j*dstack);
        /**average the rfs*/ 
            for (i=1;i<mklen/2;++i){
                sq_ampc=mrcorr[2*i]*mrcorr[2*i]+mrcorr[2*i+1]*mrcorr[2*i+1];
                //sq_ampr=mrrf[2*i]*mrrf[2*i]+mrrf[2*i+1]*mrrf[2*i+1];
                //variance=(1-sq_ampc)/(float)(mnwin-1)/sq_ampc*sq_ampr;
                variance=sq_ampc;
                rf_ave[2*i+j*mklen1] += mrrf[2*i]*variance;
                rf_ave[2*i+1+j*mklen1] += mrrf[2*i+1]*variance;
                averages[2*i+j*mklen1] += variance;
                averages[2*i+1+j*mklen1] += variance;
                }

        /* zero frequency*/
            sq_ampc=mrcorr[0]*mrcorr[0];
            //sq_ampr=mrrf[0]*mrrf[0];
            //variance=(1-sq_ampc)/(float)(mnwin-1)/sq_ampc*sq_ampr;
            variance=sq_ampc;
            rf_ave[0+j*mklen1] += mrrf[0]*variance; averages[0+j*mklen1] += variance;

        /*nyqiust frequency*/
            variance=1;
            rf_ave[1+j*mklen1] += mrrf[1]*variance; averages[1+j*mklen1] += variance;

            counts_stack[j]+=1;
        }
        fprintf(stderr,"\n");
    }
     
    for (i=0;i<nstack;++i){
        if (counts_stack[i]==0) continue;
        fprintf(stderr,"%d-----%d:%f    counts:%d\n",i,stack_t,minstack+i*dstack,counts_stack[i]);
        for (j=0;j<mklen1;++j)  rf_ave[j+i*mklen1]=rf_ave[j+i*mklen1]/averages[j+i*mklen1]; 
        sign=-1;
        jrealft(rf_ave+i*mklen1-1,(unsigned long)mklen1,sign);
        for (j=0;j<mklen1;++j) rf_ave[i*mklen1+j]=rf_ave[i*mklen1+j]*2/mklen1;
        hdz.npts=mklen1;
        hdz.user9=counts_stack[i];
        hdz.user8=minstack+i*dstack;
        sprintf(sacout,"%s/rf_stack_%d_%f.sac",outdir,stack_t,minstack+i*dstack);
        write_sac(sacout,hdz,rf_ave+i*mklen1);
    }
     
    return 0;

}
