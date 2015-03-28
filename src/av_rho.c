/* $Id: av_rho.c,v 1.1.1.1 2007/11/21 22:44:52 koehn Exp $*/

#include "fd.h"

void av_rho(float **rho, float **rip, float **rjp){

	extern int NX, NY;
	int i, j;
		
	
	for (j=1;j<=NY;j++){
		for (i=1;i<=NX;i++){
	       
			rjp[j][i] = 1.0/(0.5*(rho[j][i]+rho[j+1][i])); 	
			rip[j][i] = 1.0/(0.5*(rho[j][i]+rho[j][i+1])); 	
	
		}
	}
}
