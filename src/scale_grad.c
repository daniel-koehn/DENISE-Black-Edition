/*
 * Scale Matrix A by factor a
 *
 * Daniel Koehn
 * Kiel, 22/08/2017
 */

#include "fd.h"

void scale_grad(float ** A, float a, float ** B, int n, int m){

	/* local variables */
	int i, j;
	
	/* scale matrix A by factor a */
	for (j=1;j<=m;j++){
	   for (i=1;i<=n;i++){
	            
	       B[j][i] = a * A[j][i];
	            
	    }
	}

}



