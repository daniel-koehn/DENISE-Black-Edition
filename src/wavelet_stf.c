/*------------------------------------------------------------------------
*   Passing the soource time function of stf.c to signals
*   M. Schaefer Sept. 2011 KIT
*  ----------------------------------------------------------------------*/

#include "fd.h"


float ** wavelet_stf(int nsrc, int ishot, float ** signals_stf){


	/* extern variables */
	extern int QUELLART, NT, MYID, INV_STF;
	extern float  DT;
	extern FILE *FP;

	/*local variables */
	int nt, k, z=1;
	float ** signals;
	
	signals=fmatrix(1,nsrc,1,NT);
		
		for (k=1;k<=nsrc;k++){
			for (nt=1;nt<=NT;nt++) {
						signals[z][nt]=signals_stf[k][nt];	
						}
					++z;	
					}
	
	fprintf(FP," Message from function wavelet_stf written by PE %d \n",MYID);
	fprintf(FP," %d source positions located in subdomain of PE %d \n",nsrc,MYID);
	fprintf(FP," have been assigned with a source signal out of source time function. \n");
			
	return signals;	

}
