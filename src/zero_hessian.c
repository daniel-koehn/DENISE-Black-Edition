/*------------------------------------------------------------------------
 *   zero wavefield
 *  
 *  
 *   last update 17/02/07, T. Bohlen
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void zero_hessian(int ny1, int ny2, int nx1, int nx2, int nshots, float *** green_vx, float *** greeni_vx, float *** green_vy, float *** greeni_vy, float *** green_sxx, float *** greeni_sxx, float *** green_syy, 
                 float *** greeni_syy, float *** green_sxy, float *** greeni_sxy){

	register int i, j, k;
	extern int FW, NX, NY;

	
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=1;k<=nshots;k++){
			
				green_vx[j][i][k]=0.0;
				green_vy[j][i][k]=0.0;
				green_sxx[j][i][k]=0.0;
				green_syy[j][i][k]=0.0;
                                green_sxy[j][i][k]=0.0;
				
				greeni_vx[j][i][k]=0.0;
				greeni_vy[j][i][k]=0.0;
				greeni_sxx[j][i][k]=0.0;
				greeni_syy[j][i][k]=0.0;
				greeni_sxy[j][i][k]=0.0;
				
				}
			}
		}           
	
}
