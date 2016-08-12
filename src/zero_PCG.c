/*------------------------------------------------------------------------
 *   zero PCG-vectors
 *  
 *  
 *   last update 12/08/2016, D. Koehn
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void zero_PCG(float * PCG_old, float * PCG_new, float * PCG_dir, int PCG_vec){

	register int i;

	for (i=1;i<=PCG_vec;i++){
            PCG_old[i]=0.0;
            PCG_new[i]=0.0;
            PCG_dir[i]=0.0;
	}
	
	
}
