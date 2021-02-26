#include "jl.h"


int adwait(double *sqr_spec,  double *dcf,
            double *el, int nwin, int num_freq, double *ares, double *degf, double avar)
{
/*  input:
sqr_spec = squared modulus of eigenspectra
el = eigenvalues
nwin = number of taper windows
num_freq = number of frequency points

return:
   ares = spectrum estimate
   degf = degrees of freedom
   dcf  = 
          
  this version uses thomson's algorithm (p. 1065) for calculating 
  the adaptive spectrum estimate, 
*/
      double as,das,tol,a1,scale,ax,fn,fx;
      double *spw, *bias;
      double test_tol, dif;
      int jitter, i,j,k, kpoint, jloop;
       float df;
       /*
c  set tolerance for iterative scheme exit */


      tol=3.0e-4;
      jitter=0;
      scale=avar;
                 /***********************************
                   we scale the bias by the total variance of the frequency transform
                   from zero freq to the nyquist
                  in this application we scale the eigenspectra by the bias in order to avoid
                   possible floating point overflow
                 ************************************/
	spw = (double *) malloc((size_t) nwin * sizeof(double));
	bias = (double *) malloc((size_t) nwin * sizeof(double));

      
      for( i=0;i<nwin; i++)
          {
            
            bias[i]=(1.00-el[i]);
             }


       /* START */
    for( jloop=0; jloop<num_freq; jloop++)
    {   
        
       for( i=0;i<nwin; i++)
         {  kpoint=jloop+i*num_freq;
            spw[i]=(sqr_spec[kpoint])/scale ;
           }
                          /********************************************
                            first guess is the average of the two 
                              lowest-order eigenspectral estimates
                           ********************************************/
       as=(spw[0]+spw[1])/2.00;

                              
                              /*   find coefficients */

        for( k=0; k<20 ; k++) 
        {
          fn=0.00;
          fx=0.00;

          for( i=0;i<nwin; i++)
           {
               a1=sqrt(el[i])*as/(el[i]*as+bias[i]);
               a1=a1*a1;
               fn=fn+a1*spw[i];
               fx=fx+a1;
           }
  

         ax=fn/fx;
         dif = ax-as;
         das=ABS(dif);
      
         test_tol = das/as;
         if( test_tol < tol )
            { 
                   break;
               }

         as=ax;
        }

 
                           /*   flag if iteration does not converge */

      if(k>=20)  jitter++;
  
       ares[jloop]=as*scale;
                            /*    calculate degrees of freedom */
      df=0.0;
      for( i=0;i< nwin; i++)
       {
          kpoint=jloop+i*num_freq;
          dcf[kpoint]=sqrt(el[i])*as/(el[i]*as+bias[i]);
          df=df+dcf[kpoint]*dcf[kpoint];
       }
 			/*
			 * we normalize degrees of freedom by the weight of
			 * the first eigenspectrum this way we never have
			 * fewer than two degrees of freedom
			 */

       degf[jloop]=df*2./(dcf[jloop]*dcf[jloop]);

  }                                       /* END */

     fprintf(stderr,"%d failed iterations\n",jitter);
      free(spw);
      free(bias);

     return jitter;
}
