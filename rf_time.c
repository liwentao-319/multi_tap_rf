/**********************************************************************************
 *    this program is to calculate receiver function time series using data from 
 * program rf_spec. 
 * 
 * 
 * 
 * *******************************************************************************/
#include "jl.h"
#include "sac.h"

void jrealft(float data[], unsigned long n, int isign);
int is_pow_2(unsigned long iii);
void savedata(char *filename,float *data,unsigned long npts);
int main(int argc, char **argv){

    char filein[128],fileout[128],flag;
    int i,j,isiof=0, sign=-1,pow2=1,error=0;
    unsigned long klen;
    float *data;
    float fc,df,ff,dt,cosin;
    float *rff_raw,*rff_ed;////
    SACHEAD hdi;
    FILE *fd;
    fc=2;
    dt=0.01;
    /*input argument*/
    for (i=1; !error && i < argc; i++) {
        if (argv[i][0] == '-') {
            switch(argv[i][1]) {

            case 'I':  
                sscanf(&argv[i][2],"%s",filein);
                isiof += 1;
                break;
    
            case 'O':
                sscanf(&argv[i][2],"%s",fileout);
                isiof += 5;
                break;
    
            case 'F':
                sscanf(&argv[i][2],"%f",&fc);
                isiof += 7;
                break;

            default:
                error = TRUE;
                break;
        }
        } 
    }

    if (isiof!=13) {
        fprintf(stderr,"error message:%i\n",isiof);
        fprintf(stderr,"usage: %s -Isrc -Oout -Ffc\n\
        -I  filename of RF spectra to load\n\
        -O  filename of RF time series to save\n\
        -F  corner frequency of lowpass filter with cos^2(pi*f/fc/2)\n",argv[0]);
        return -1;
    }
    
    data = read_sac(filein,&hdi);
    klen = hdi.npts;
    dt = hdi.delta;
    fprintf(stderr,"filein:%s\n\
     fileout:%s\n\
     fc=%f,dt=%f,klen=%ld\n",filein,fileout,fc,dt,klen);
    pow2=is_pow_2(klen);
    if (pow2==0){
        fprintf(stderr,"wrong number of points input(%ld), which must be pow2 ",klen);
        return -1;
    }

    /* add low_pass filter to rf with cos^2(pi*f/fc/2)*/
    /*save spectra of rff with filtered and without filtered*/
rff_raw=(float *)malloc((size_t)klen/2*sizeof(float));///
rff_ed=(float *)malloc((size_t)klen/2*sizeof(float));///

    df=1/dt/klen;
    for (i=1;i<klen/2;++i){
        ff=df*i;
        if (ff<=fc)  cosin=(float)cos((double) PI*ff/fc/2);
        else cosin=0;
rff_raw[i]=data[2*i]*data[2*i]+data[2*i+1]*data[2*i+1];/////
        data[2*i]=data[2*i]*cosin*cosin;
        data[2*i+1]=data[2*i+1]*cosin*cosin;
rff_ed[i]=data[2*i]*data[2*i]+data[2*i+1]*data[2*i+1];/////
    }
    ff=df*klen/2;

rff_raw[0]=data[0]*data[0];////
rff_ed[0]=data[0]*data[0];////
rff_raw[1]=data[1]*data[1];////

    /**freq=nyquist frequecy*/
    if (ff<=fc)   cosin=(float)cos((double) PI*ff/fc/2);
    else  cosin=0;
    data[1]=data[1]*cosin*cosin;

rff_ed[1]=data[1] ;////
savedata("rff_raw.txt",rff_raw,(unsigned long)klen/2);
savedata("rff_ed.txt",rff_ed,(unsigned long)klen/2);

    /* reverse fft*/
    jrealft(data-1,klen,sign);
    /*save data*/
    write_sac(fileout,hdi,data);


    

}
