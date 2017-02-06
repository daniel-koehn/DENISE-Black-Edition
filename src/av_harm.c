/*-------------------------------------------
 *  Harmonic averaging of model parameter m  
 *
 *  D. Koehn
 *  Kiel, 02.02.2017
 *  ----------------------------------------- */
#include "fd.h"

void av_harm(float ** m, float ** mh){

	extern int NX, NY;
	int i, j;	
	
	for (j=1;j<=NY;j++){
		for (i=1;i<=NX;i++){
	       
		       mh[j][i]=4.0/((1.0/m[j][i])+(1.0/m[j][i+1])+(1.0/m[j+1][i])+(1.0/m[j+1][i+1])); 
				
		       if((m[j][i]==0.0)||(m[j][i+1]==0.0)||(m[j+1][i]==0.0)||(m[j+1][i+1]==0.0)){
		           mh[j][i]=0.0;
		       }
		      		
		}
	}	
	
}
