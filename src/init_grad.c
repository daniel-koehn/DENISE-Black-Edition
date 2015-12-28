/*
 * Initialize gradient 
 *
 * Daniel Koehn
 * Kiel, 11/12/2015
 */

#include "fd.h"

void init_grad(float ** A){

        /* global variables */
	extern int NX, NY, INVMAT;

	/* local variables */
	int i, j;
	
	/* initiate gradient A */
	for (i=1;i<=NX;i++){
	    for (j=1;j<=NY;j++){
		    
	       A[j][i] = 0.0;
		    
	    }
	}

}



