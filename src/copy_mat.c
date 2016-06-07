/*
 * Copy matrix A to matrix B
 *
 * Daniel Koehn
 * Kiel, 03/06/2016
 */

#include "fd.h"

void copy_mat(float ** A, float ** B){

        /* global variables */
	extern int NX, NY;

	/* local variables */
	int i, j;
	
	/* copy matrix A to matrix B */
	for (i=1;i<=NX;i++){
	    for (j=1;j<=NY;j++){
		    
	       B[j][i] = A[j][i];
		    
	    }
	}

}



