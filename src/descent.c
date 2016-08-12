/*------------------------------------------------------------------------
 * Calculate steepest descent direction from gradient
 *
 * 
 * Daniel Koehn
 * Kiel, 14.12.2015
 * ----------------------------------------------------------------------*/

#include "fd.h"

void descent(float ** grad, float ** gradm){

        /* global variables */
	extern int NX, NY;
	
        /* local variables */
	int i, j;

        for (i=1;i<=NX;i++){
             for (j=1;j<=NY;j++){
   	  
	          gradm[j][i] = -grad[j][i];

             }
        }

}



