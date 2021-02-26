#include "jl.h"
#include<sys/stat.h>
#include<sys/types.h>
#include <unistd.h>


/***************************************************************************/
int
get_pow_2(int inum)
{
	int             j, klength;
	/* find smallest power of 2 that encompasses the data */

	for (j = 1; pow((double) 2, (double) j) < inum; j++);
	return klength = pow((double) 2, (double) j);
}


double  remove_mean(float x[], int lx)
   {

     int k;
      double mean;
      mean = 0.;
     if(lx < 2 ) return mean;
 
      for( k=0; k<=lx; k++)
        {
        mean = x[k]+mean;
        }

       mean = mean/ (float)lx;

      for( k=0; k<=lx; k++)
        {
        x[k] = x[k] - mean;
        }

      
      return  mean;
   }
/*********************************************************/
void zero_pad(float  output[], int start , int olength)
{
    int i;
    for( i= start ; i< olength; i++) 
         {   
            output[i] = 0.0; 
               }
}
/***********************************************************************/
float get_cos_taper(int n, int k)
{
int l;
float vwin;
vwin = 0.0;

if(k<0 || k > n) return vwin;
vwin = 1.0;

      l=(n-2)/10;
      if(k<=l) vwin=0.5*(1.0-cos(k*PI/(l+1)));
      if(k>=n-l-2) vwin=0.5*(1.0-cos((n-k-1)*PI/(l+1)));


        return vwin;
}
/**********************************************************/
int is_pow_2(unsigned long iii){
  if (iii<=1) return 0;
  while(iii>1){
    if(iii%2!=0) return 0;
    iii=iii/2;
  }
  return 1;


}

/***************************************************************/

void savedata(char *filename,float *data,unsigned long npts){
    FILE *fd;
    int i;
    if ((fd=fopen(filename,"w"))==NULL){
    printf("error when created a new file\n");
    }
    for (i=0;i<npts;++i){
        fprintf(fd,"%f\n",data[i]);

    }
    fclose(fd);
}
/*****************************************************************/
int strrep(char * str1, char * str2, char * str3){
    int i, j, k, done, count = 0, gap = 0;
    char temp[300];
    for(i = 0; i < strlen(str1); i += gap){
        if(str1[i] == str2[0]){
            done = 0;
            for(j = i, k = 0; k < strlen(str2); j++, k++){
                if(str1[j] != str2[k]){
                    done = 1;
                    gap = k;
                    break;
                }
            }
            if(done == 0){ // 已找到待替换字符串并替换
                for(j = i + strlen(str2), k = 0; j < strlen(str1); j++, k++){ // 保存原字符串中剩余的字符
                    temp[k] = str1[j];
                }
                temp[k] = '\0'; // 将字符数组变成字符串
                for(j = i, k = 0; k < strlen(str3); j++, k++){ // 字符串替换
                    str1[j] = str3[k];
                    count++;
                }
                for(k = 0; k < strlen(temp); j++, k++){ // 剩余字符串回接
                    str1[j] = temp[k];
                }
                str1[j] = '\0'; // 将字符数组变成字符串
                gap = strlen(str2);
            }
        }else{
            gap = 1;
        }
    }
    return count;
    
}


/********************************************************/
void cos2_oneside_win(float dt, float fc, float *data, int klen){
    float cosin,ff,df=1/dt/klen;
    int ii;
    for (ii=1;ii<klen/2;++ii){
        ff=df*ii;
        if (ff<=fc)  cosin=(float)cos((double) PI*ff/fc/2);
        else cosin=0;
        data[2*ii]=data[2*ii]*cosin*cosin;
        data[2*ii+1]=data[2*ii+1]*cosin*cosin;
    }
    ff=df*klen/2;
    /**freq=nyquist frequecy*/
    if (ff<=fc)   cosin=(float)cos((double) PI*ff/fc/2);
    else  cosin=0;
    data[1]=data[1]*cosin*cosin;
}


/*****************************************************************/

void mkdirs(char *muldir) 
{
    int i,len;
    char str[512];    
    strncpy(str, muldir, 512);
    len=strlen(str);
    for( i=0; i<len; i++ )
    {
        if( str[i]=='/' )
        {
            str[i] = '\0';
            if( access(str,0)!=0 )
            {
                mkdir( str, 0777 );
            }
            str[i]='/';
        }
    }
    if( len>0 && access(str,0)!=0 )
    {
        mkdir( str, 0777 );
    }
    return;
}