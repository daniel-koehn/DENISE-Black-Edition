/*------------------------------------------------------------------------
 *   generate P-wave source at source nodes
 *   last update 12/02/02, T. Bohlen
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void psource(int nt, float ** sxx, float ** syy,
float **  srcpos_loc, float ** signals, int nsrc){


        extern int MYID, NT;
	extern float DH, DT;
	int i, j, l;
	float amp;



	/* adding source wavelet to stress components 
	   (explosive source) at source points */


	for (l=1;l<=nsrc;l++) {
		i=(int)srcpos_loc[1][l];
		j=(int)srcpos_loc[2][l];
		
		if(nt==1){amp=signals[l][nt+1]/DT;}
		if((nt>1)&&(nt<NT)){amp=(signals[l][nt+1]-signals[l][nt-1])/DT;}
		if(nt==NT){amp=-signals[l][nt-1]/DT;}
		
		/*amp = signals[l][nt];*/
		
		sxx[j][i]+=amp;
		syy[j][i]+=amp;
	}

}
