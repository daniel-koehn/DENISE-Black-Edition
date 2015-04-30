/*------------------------------------------------------------------------
 *   Source defined by moment tensor
 *   Daniel Koehn
 *   Kiel, the 30th of April 2015
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void msource(int nt, float ** sxx, float ** syy, float ** sxy, float **  srcpos_loc, float ** signals, int nsrc, int sw){


        extern int MYID, NT;
	extern float DH, DT;
	int i, j, l;
	float amp, M11, M12, M22;

        M11 = 0.0;
        M12 = 1.0;
        M22 = 0.0;


	/* adding source wavelet to stress components (moment-tensor source) at source points according to Gassner (2014) */

	   for (l=1;l<=nsrc;l++) {

	   	   i=(int)srcpos_loc[1][l];
		   j=(int)srcpos_loc[2][l];
		
		   if(nt==1){amp=signals[l][nt+1]/DT;}
		   if((nt>1)&&(nt<NT)){amp=(signals[l][nt+1]-signals[l][nt-1])/DT;}
		   if(nt==NT){amp=-signals[l][nt-1]/DT;}
		
		   sxx[j][i]+=amp*M11;
		   syy[j][i]+=amp*M22;
                   sxy[j][i]+=0.25*amp*M12;
                   sxy[j][i-1]+=0.25*amp*M12;
                   sxy[j-1][i]+=0.25*amp*M12;
                   sxy[j-1][i-1]+=0.25*amp*M12;
                   
	   }

}
