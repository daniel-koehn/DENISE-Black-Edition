/*-----------------------------
 *  Harmonic averaging of c44  
 *
 *  D. Koehn
 *  Kiel, 02.02.2017
 *  --------------------------- */
#include "fd.h"

void av_c44(float ** c44, float ** c44h){

	extern int NX, NY;
	int i, j;	
	
	for (j=1;j<=NY;j++){
		for (i=1;i<=NX;i++){
	       
		       c44h[j][i]=4.0/((1.0/c44[j][i])+(1.0/c44[j][i+1])+(1.0/c44[j+1][i])+(1.0/c44[j+1][i+1])); 
				
		       if((c44[j][i]==0.0)||(c44[j][i+1]==0.0)||(c44[j+1][i]==0.0)||(c44[j+1][i+1]==0.0)){
		           c44h[j][i]=0.0;
		       }
		      	
	
		}
	}	
	
}
