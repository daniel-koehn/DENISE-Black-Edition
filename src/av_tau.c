/* $Id: av_tau.c,v 2.4 2007/08/21 13:16:19 tbohlen Exp $*/

#include "fd.h"

void av_tau(float **taus, float **tausipjp){

	extern int NX, NY;
	int i, j;
	for (j=1;j<=NY;j++){
		for (i=1;i<=NX;i++){

		       tausipjp[j][i] = 0.25*(taus[j][i]+taus[j][i+1]+taus[j+1][i]+taus[j+1][i+1]);
		}
	}
}
