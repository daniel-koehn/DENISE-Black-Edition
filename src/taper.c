/*------------------------------------------------------------------------
 *   Write seismograms to disk                                  
 *   last update 19/01/02, T. Bohlen
 *  ----------------------------------------------------------------------*/
#include "fd.h"
#include "segy.h"

void  taper(float **sectionpdiff, int ntr, int ns){

	/* declaration of extern variables */
	extern int NDT, MYID, TAPERLENGTH;
	extern char DATA_DIR[STRING_SIZE];
	extern float  TIME, DH, DT, REFREC[4];
	extern int TAPERLENGTH;
        char data[STRING_SIZE];
	const float xshift=800.0, yshift=800.0;
	
	/* declaration of local variables */
	int i,j, h;
	segy tr;
	int tracl1 ;
	float xr, yr, x, y, dump, a;
	float damping, amp;
	float *window, *amp1;
	
	window = vector(1,ntr);
        amp1 = vector(1,ntr);
        
	/* define taper function */
	/* --------------------- */

	/* Hanning window */  

        /*for (i=1;i<=ntr;i++){
	window[i] = 1.0;
	
	    if(i<=TAPERLENGTH/2){
	       window[i] = 0.5 + 0.5 * cos(2*PI*(i+(TAPERLENGTH/2))/TAPERLENGTH);
	    }
	    
	    if(i>=(ntr-(TAPERLENGTH/2))){
	       window[i] = 0.5 + 0.5 * cos(2*PI*(i-(ntr+TAPERLENGTH/2))/TAPERLENGTH);
	    }
	}*/

	/* Hamming window */  

        /*for (i=1;i<=ntr;i++){
	window[i] = 1.0;
	
	    if(i<=TAPERLENGTH/2){
	       window[i] = 0.54 + 0.46 * cos(2*PI*(i+(TAPERLENGTH/2))/TAPERLENGTH);
	    }
	    
	    if(i>=(ntr-(TAPERLENGTH/2))){
	       window[i] = 0.54 + 0.46 * cos(2*PI*(i-(ntr+TAPERLENGTH/2))/TAPERLENGTH);
	    }
	}*/
	
	/* Gaussian window */  
        a=3;    /* damping coefficient */
        /*for (i=1;i<=ntr;i++){
	window[i] = 1.0;
	
	    if(i<=TAPERLENGTH){
	       window[i] = exp(-(1.0/2.0)*(a*(i-TAPERLENGTH)/(TAPERLENGTH/2.0))*(a*(i-TAPERLENGTH)/(TAPERLENGTH/2.0)));
	    }
	    
	    if(i>=(ntr-TAPERLENGTH)){
	       window[i] = exp(-(1.0/2.0)*(a*(i-(ntr-TAPERLENGTH))/(TAPERLENGTH/2.0))*(a*(i-(ntr-TAPERLENGTH))/(TAPERLENGTH/2.0)));
	    }
	    
	}*/
	
	/* Poisson window */  
        a=3;    /* damping coefficient */
        /*for (i=1;i<=ntr;i++){
	window[i] = 1.0;
	
	    if(i<=TAPERLENGTH){
	       window[i] = exp(-a*abs(i-TAPERLENGTH)/(TAPERLENGTH/2.0));
	    }
	    
	    if(i>=(ntr-TAPERLENGTH)){
	       window[i] = exp(-a*abs(i-(ntr-TAPERLENGTH))/(TAPERLENGTH/2.0));
	    }
	    
	}*/
	
	/* "Cerjan"-Window */
        damping=90.0;
        amp=1.0-damping/100.0;
        a=sqrt(-log(amp)/((TAPERLENGTH-1)*(TAPERLENGTH-1)));
        for (i=1;i<=ntr;i++){window[i]=1.0;}
        for (i=1;i<=TAPERLENGTH;i++){amp1[i]=exp(-(a*a*(TAPERLENGTH-i)*(TAPERLENGTH-i)));}
        for (i=1;i<=TAPERLENGTH;i++){window[i]=amp1[i];}
        h=1;
        for (i=ntr;i>=(ntr-TAPERLENGTH+3);i--){window[i]=amp1[h];h++;}
	
	/* Apply taper window */
		for(tracl1=1;tracl1<=ntr;tracl1++){ 
                        
			for(j=1;j<=ns;j++){
			  sectionpdiff[tracl1][j+1]*=window[tracl1];
			}
		}

}
